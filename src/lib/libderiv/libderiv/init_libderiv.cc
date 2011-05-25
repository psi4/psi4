#include <stdlib.h>
#include <strings.h>
#include <libint/libint.h>
#include "libderiv.h"
#include "d1hrr_header.h"

extern "C" {
void (*build_deriv1_eri[5][5][5][5])(Libderiv_t *, int);

void (*build_deriv12_eri[4][4][4][4])(Libderiv_t *, int);

int libderiv1_stack_size[5];
int libderiv12_stack_size[4];
void init_libderiv_base()
{
  build_deriv1_eri[0][0][0][0] = d1hrr_order_0000;
  build_deriv1_eri[0][0][1][0] = d1hrr_order_00p0;
  build_deriv1_eri[0][0][2][0] = d1hrr_order_00d0;
  build_deriv1_eri[0][0][1][1] = d1hrr_order_00pp;
  build_deriv1_eri[0][0][3][0] = d1hrr_order_00f0;
  build_deriv1_eri[0][0][2][1] = d1hrr_order_00dp;
  build_deriv1_eri[0][0][4][0] = d1hrr_order_00g0;
  build_deriv1_eri[0][0][3][1] = d1hrr_order_00fp;
  build_deriv1_eri[0][0][2][2] = d1hrr_order_00dd;
  build_deriv1_eri[0][0][4][1] = d1hrr_order_00gp;
  build_deriv1_eri[0][0][3][2] = d1hrr_order_00fd;
  build_deriv1_eri[0][0][4][2] = d1hrr_order_00gd;
  build_deriv1_eri[0][0][3][3] = d1hrr_order_00ff;
  build_deriv1_eri[0][0][4][3] = d1hrr_order_00gf;
  build_deriv1_eri[0][0][4][4] = d1hrr_order_00gg;
  build_deriv1_eri[1][0][1][0] = d1hrr_order_p0p0;
  build_deriv1_eri[1][0][2][0] = d1hrr_order_p0d0;
  build_deriv1_eri[1][0][1][1] = d1hrr_order_p0pp;
  build_deriv1_eri[1][0][3][0] = d1hrr_order_p0f0;
  build_deriv1_eri[1][0][2][1] = d1hrr_order_p0dp;
  build_deriv1_eri[1][0][4][0] = d1hrr_order_p0g0;
  build_deriv1_eri[1][0][3][1] = d1hrr_order_p0fp;
  build_deriv1_eri[1][0][2][2] = d1hrr_order_p0dd;
  build_deriv1_eri[1][0][4][1] = d1hrr_order_p0gp;
  build_deriv1_eri[1][0][3][2] = d1hrr_order_p0fd;
  build_deriv1_eri[1][0][4][2] = d1hrr_order_p0gd;
  build_deriv1_eri[1][0][3][3] = d1hrr_order_p0ff;
  build_deriv1_eri[1][0][4][3] = d1hrr_order_p0gf;
  build_deriv1_eri[1][0][4][4] = d1hrr_order_p0gg;
  build_deriv1_eri[2][0][2][0] = d1hrr_order_d0d0;
  build_deriv1_eri[2][0][1][1] = d1hrr_order_d0pp;
  build_deriv1_eri[2][0][3][0] = d1hrr_order_d0f0;
  build_deriv1_eri[2][0][2][1] = d1hrr_order_d0dp;
  build_deriv1_eri[2][0][4][0] = d1hrr_order_d0g0;
  build_deriv1_eri[2][0][3][1] = d1hrr_order_d0fp;
  build_deriv1_eri[2][0][2][2] = d1hrr_order_d0dd;
  build_deriv1_eri[2][0][4][1] = d1hrr_order_d0gp;
  build_deriv1_eri[2][0][3][2] = d1hrr_order_d0fd;
  build_deriv1_eri[2][0][4][2] = d1hrr_order_d0gd;
  build_deriv1_eri[2][0][3][3] = d1hrr_order_d0ff;
  build_deriv1_eri[2][0][4][3] = d1hrr_order_d0gf;
  build_deriv1_eri[2][0][4][4] = d1hrr_order_d0gg;
  build_deriv1_eri[1][1][2][0] = d1hrr_order_ppd0;
  build_deriv1_eri[1][1][1][1] = d1hrr_order_pppp;
  build_deriv1_eri[1][1][3][0] = d1hrr_order_ppf0;
  build_deriv1_eri[1][1][2][1] = d1hrr_order_ppdp;
  build_deriv1_eri[1][1][4][0] = d1hrr_order_ppg0;
  build_deriv1_eri[1][1][3][1] = d1hrr_order_ppfp;
  build_deriv1_eri[1][1][2][2] = d1hrr_order_ppdd;
  build_deriv1_eri[1][1][4][1] = d1hrr_order_ppgp;
  build_deriv1_eri[1][1][3][2] = d1hrr_order_ppfd;
  build_deriv1_eri[1][1][4][2] = d1hrr_order_ppgd;
  build_deriv1_eri[1][1][3][3] = d1hrr_order_ppff;
  build_deriv1_eri[1][1][4][3] = d1hrr_order_ppgf;
  build_deriv1_eri[1][1][4][4] = d1hrr_order_ppgg;
  build_deriv1_eri[3][0][3][0] = d1hrr_order_f0f0;
  build_deriv1_eri[3][0][2][1] = d1hrr_order_f0dp;
  build_deriv1_eri[3][0][4][0] = d1hrr_order_f0g0;
  build_deriv1_eri[3][0][3][1] = d1hrr_order_f0fp;
  build_deriv1_eri[3][0][2][2] = d1hrr_order_f0dd;
  build_deriv1_eri[3][0][4][1] = d1hrr_order_f0gp;
  build_deriv1_eri[3][0][3][2] = d1hrr_order_f0fd;
  build_deriv1_eri[3][0][4][2] = d1hrr_order_f0gd;
  build_deriv1_eri[3][0][3][3] = d1hrr_order_f0ff;
  build_deriv1_eri[3][0][4][3] = d1hrr_order_f0gf;
  build_deriv1_eri[3][0][4][4] = d1hrr_order_f0gg;
  build_deriv1_eri[2][1][3][0] = d1hrr_order_dpf0;
  build_deriv1_eri[2][1][2][1] = d1hrr_order_dpdp;
  build_deriv1_eri[2][1][4][0] = d1hrr_order_dpg0;
  build_deriv1_eri[2][1][3][1] = d1hrr_order_dpfp;
  build_deriv1_eri[2][1][2][2] = d1hrr_order_dpdd;
  build_deriv1_eri[2][1][4][1] = d1hrr_order_dpgp;
  build_deriv1_eri[2][1][3][2] = d1hrr_order_dpfd;
  build_deriv1_eri[2][1][4][2] = d1hrr_order_dpgd;
  build_deriv1_eri[2][1][3][3] = d1hrr_order_dpff;
  build_deriv1_eri[2][1][4][3] = d1hrr_order_dpgf;
  build_deriv1_eri[2][1][4][4] = d1hrr_order_dpgg;
  build_deriv1_eri[4][0][4][0] = d1hrr_order_g0g0;
  build_deriv1_eri[4][0][3][1] = d1hrr_order_g0fp;
  build_deriv1_eri[4][0][2][2] = d1hrr_order_g0dd;
  build_deriv1_eri[4][0][4][1] = d1hrr_order_g0gp;
  build_deriv1_eri[4][0][3][2] = d1hrr_order_g0fd;
  build_deriv1_eri[4][0][4][2] = d1hrr_order_g0gd;
  build_deriv1_eri[4][0][3][3] = d1hrr_order_g0ff;
  build_deriv1_eri[4][0][4][3] = d1hrr_order_g0gf;
  build_deriv1_eri[4][0][4][4] = d1hrr_order_g0gg;
  build_deriv1_eri[3][1][4][0] = d1hrr_order_fpg0;
  build_deriv1_eri[3][1][3][1] = d1hrr_order_fpfp;
  build_deriv1_eri[3][1][2][2] = d1hrr_order_fpdd;
  build_deriv1_eri[3][1][4][1] = d1hrr_order_fpgp;
  build_deriv1_eri[3][1][3][2] = d1hrr_order_fpfd;
  build_deriv1_eri[3][1][4][2] = d1hrr_order_fpgd;
  build_deriv1_eri[3][1][3][3] = d1hrr_order_fpff;
  build_deriv1_eri[3][1][4][3] = d1hrr_order_fpgf;
  build_deriv1_eri[3][1][4][4] = d1hrr_order_fpgg;
  build_deriv1_eri[2][2][4][0] = d1hrr_order_ddg0;
  build_deriv1_eri[2][2][3][1] = d1hrr_order_ddfp;
  build_deriv1_eri[2][2][2][2] = d1hrr_order_dddd;
  build_deriv1_eri[2][2][4][1] = d1hrr_order_ddgp;
  build_deriv1_eri[2][2][3][2] = d1hrr_order_ddfd;
  build_deriv1_eri[2][2][4][2] = d1hrr_order_ddgd;
  build_deriv1_eri[2][2][3][3] = d1hrr_order_ddff;
  build_deriv1_eri[2][2][4][3] = d1hrr_order_ddgf;
  build_deriv1_eri[2][2][4][4] = d1hrr_order_ddgg;
  build_deriv1_eri[4][1][4][1] = d1hrr_order_gpgp;
  build_deriv1_eri[4][1][3][2] = d1hrr_order_gpfd;
  build_deriv1_eri[4][1][4][2] = d1hrr_order_gpgd;
  build_deriv1_eri[4][1][3][3] = d1hrr_order_gpff;
  build_deriv1_eri[4][1][4][3] = d1hrr_order_gpgf;
  build_deriv1_eri[4][1][4][4] = d1hrr_order_gpgg;
  build_deriv1_eri[3][2][4][1] = d1hrr_order_fdgp;
  build_deriv1_eri[3][2][3][2] = d1hrr_order_fdfd;
  build_deriv1_eri[3][2][4][2] = d1hrr_order_fdgd;
  build_deriv1_eri[3][2][3][3] = d1hrr_order_fdff;
  build_deriv1_eri[3][2][4][3] = d1hrr_order_fdgf;
  build_deriv1_eri[3][2][4][4] = d1hrr_order_fdgg;
  build_deriv1_eri[4][2][4][2] = d1hrr_order_gdgd;
  build_deriv1_eri[4][2][3][3] = d1hrr_order_gdff;
  build_deriv1_eri[4][2][4][3] = d1hrr_order_gdgf;
  build_deriv1_eri[4][2][4][4] = d1hrr_order_gdgg;
  build_deriv1_eri[3][3][4][2] = d1hrr_order_ffgd;
  build_deriv1_eri[3][3][3][3] = d1hrr_order_ffff;
  build_deriv1_eri[3][3][4][3] = d1hrr_order_ffgf;
  build_deriv1_eri[3][3][4][4] = d1hrr_order_ffgg;
  build_deriv1_eri[4][3][4][3] = d1hrr_order_gfgf;
  build_deriv1_eri[4][3][4][4] = d1hrr_order_gfgg;
  build_deriv1_eri[4][4][4][4] = d1hrr_order_gggg;
  build_deriv12_eri[0][0][0][0] = d12hrr_order_0000;
  build_deriv12_eri[0][0][1][0] = d12hrr_order_00p0;
  build_deriv12_eri[0][0][2][0] = d12hrr_order_00d0;
  build_deriv12_eri[0][0][1][1] = d12hrr_order_00pp;
  build_deriv12_eri[0][0][3][0] = d12hrr_order_00f0;
  build_deriv12_eri[0][0][2][1] = d12hrr_order_00dp;
  build_deriv12_eri[0][0][3][1] = d12hrr_order_00fp;
  build_deriv12_eri[0][0][2][2] = d12hrr_order_00dd;
  build_deriv12_eri[0][0][3][2] = d12hrr_order_00fd;
  build_deriv12_eri[0][0][3][3] = d12hrr_order_00ff;
  build_deriv12_eri[1][0][1][0] = d12hrr_order_p0p0;
  build_deriv12_eri[1][0][2][0] = d12hrr_order_p0d0;
  build_deriv12_eri[1][0][1][1] = d12hrr_order_p0pp;
  build_deriv12_eri[1][0][3][0] = d12hrr_order_p0f0;
  build_deriv12_eri[1][0][2][1] = d12hrr_order_p0dp;
  build_deriv12_eri[1][0][3][1] = d12hrr_order_p0fp;
  build_deriv12_eri[1][0][2][2] = d12hrr_order_p0dd;
  build_deriv12_eri[1][0][3][2] = d12hrr_order_p0fd;
  build_deriv12_eri[1][0][3][3] = d12hrr_order_p0ff;
  build_deriv12_eri[2][0][2][0] = d12hrr_order_d0d0;
  build_deriv12_eri[2][0][1][1] = d12hrr_order_d0pp;
  build_deriv12_eri[2][0][3][0] = d12hrr_order_d0f0;
  build_deriv12_eri[2][0][2][1] = d12hrr_order_d0dp;
  build_deriv12_eri[2][0][3][1] = d12hrr_order_d0fp;
  build_deriv12_eri[2][0][2][2] = d12hrr_order_d0dd;
  build_deriv12_eri[2][0][3][2] = d12hrr_order_d0fd;
  build_deriv12_eri[2][0][3][3] = d12hrr_order_d0ff;
  build_deriv12_eri[1][1][2][0] = d12hrr_order_ppd0;
  build_deriv12_eri[1][1][1][1] = d12hrr_order_pppp;
  build_deriv12_eri[1][1][3][0] = d12hrr_order_ppf0;
  build_deriv12_eri[1][1][2][1] = d12hrr_order_ppdp;
  build_deriv12_eri[1][1][3][1] = d12hrr_order_ppfp;
  build_deriv12_eri[1][1][2][2] = d12hrr_order_ppdd;
  build_deriv12_eri[1][1][3][2] = d12hrr_order_ppfd;
  build_deriv12_eri[1][1][3][3] = d12hrr_order_ppff;
  build_deriv12_eri[3][0][3][0] = d12hrr_order_f0f0;
  build_deriv12_eri[3][0][2][1] = d12hrr_order_f0dp;
  build_deriv12_eri[3][0][3][1] = d12hrr_order_f0fp;
  build_deriv12_eri[3][0][2][2] = d12hrr_order_f0dd;
  build_deriv12_eri[3][0][3][2] = d12hrr_order_f0fd;
  build_deriv12_eri[3][0][3][3] = d12hrr_order_f0ff;
  build_deriv12_eri[2][1][3][0] = d12hrr_order_dpf0;
  build_deriv12_eri[2][1][2][1] = d12hrr_order_dpdp;
  build_deriv12_eri[2][1][3][1] = d12hrr_order_dpfp;
  build_deriv12_eri[2][1][2][2] = d12hrr_order_dpdd;
  build_deriv12_eri[2][1][3][2] = d12hrr_order_dpfd;
  build_deriv12_eri[2][1][3][3] = d12hrr_order_dpff;
  build_deriv12_eri[3][1][3][1] = d12hrr_order_fpfp;
  build_deriv12_eri[3][1][2][2] = d12hrr_order_fpdd;
  build_deriv12_eri[3][1][3][2] = d12hrr_order_fpfd;
  build_deriv12_eri[3][1][3][3] = d12hrr_order_fpff;
  build_deriv12_eri[2][2][3][1] = d12hrr_order_ddfp;
  build_deriv12_eri[2][2][2][2] = d12hrr_order_dddd;
  build_deriv12_eri[2][2][3][2] = d12hrr_order_ddfd;
  build_deriv12_eri[2][2][3][3] = d12hrr_order_ddff;
  build_deriv12_eri[3][2][3][2] = d12hrr_order_fdfd;
  build_deriv12_eri[3][2][3][3] = d12hrr_order_fdff;
  build_deriv12_eri[3][3][3][3] = d12hrr_order_ffff;

  libderiv1_stack_size[0] = 21;
  libderiv1_stack_size[1] = 1839;
  libderiv1_stack_size[2] = 30431;
  libderiv1_stack_size[3] = 244213;
  libderiv1_stack_size[4] = 1234656;
  libderiv12_stack_size[0] = 127;
  libderiv12_stack_size[1] = 9998;
  libderiv12_stack_size[2] = 158107;
  libderiv12_stack_size[3] = 1293745;
}

/* These functions initialize library objects */
/* Library objects operate independently of each other */
int init_libderiv1(Libderiv_t *libderiv, int max_am, int max_num_prim_quartets, int max_cart_class_size)
{
  int memory = 0;

  if (max_am >= LIBDERIV_MAX_AM1) return -1;
  libderiv->int_stack = (double *) malloc(libderiv1_stack_size[max_am]*sizeof(double));
  memory += libderiv1_stack_size[max_am];
  libderiv->zero_stack = (double *) malloc(max_cart_class_size*sizeof(double));
  bzero((char *)libderiv->zero_stack,max_cart_class_size*sizeof(double));
  memory += max_cart_class_size;
  libderiv->PrimQuartet = (prim_data *) malloc(max_num_prim_quartets*sizeof(prim_data));
  memory += max_num_prim_quartets*sizeof(prim_data)/sizeof(double);
  return memory;
}

int init_libderiv12(Libderiv_t *libderiv, int max_am, int max_num_prim_quartets, int max_cart_class_size)
{
  int memory = 0;

  if (max_am >= LIBDERIV_MAX_AM12) return -1;
  libderiv->int_stack = (double *) malloc(libderiv12_stack_size[max_am]*sizeof(double));
  memory += libderiv12_stack_size[max_am];
  libderiv->zero_stack = (double *) malloc(max_cart_class_size*sizeof(double));
  bzero((char *)libderiv->zero_stack,max_cart_class_size*sizeof(double));
  memory += max_cart_class_size;
  libderiv->PrimQuartet = (prim_data *) malloc(max_num_prim_quartets*sizeof(prim_data));
  memory += max_num_prim_quartets*sizeof(prim_data)/sizeof(double);
  return memory;
}

void free_libderiv(Libderiv_t *libderiv)
{
  if (libderiv->int_stack != NULL) {
    free(libderiv->int_stack);
    libderiv->int_stack = NULL;
  }
  if (libderiv->zero_stack != NULL) {
    free(libderiv->zero_stack);
    libderiv->zero_stack = NULL;
  }
  if (libderiv->PrimQuartet != NULL) {
    free(libderiv->PrimQuartet);
    libderiv->PrimQuartet = NULL;
  }

  return;
}

int libderiv1_storage_required(int max_am, int max_num_prim_quartets, int max_cart_class_size)
{
  int memory = 0;

  if (max_am >= LIBDERIV_MAX_AM1) return -1;
  memory += libderiv1_stack_size[max_am];
  memory += max_cart_class_size;
  memory += max_num_prim_quartets*sizeof(prim_data)/sizeof(double);
  return memory;
}
int libderiv12_storage_required(int max_am, int max_num_prim_quartets, int max_cart_class_size)
{
  int memory = 0;

  if (max_am >= LIBDERIV_MAX_AM12) return -1;
  memory += libderiv12_stack_size[max_am];
  memory += max_cart_class_size;
  memory += max_num_prim_quartets*sizeof(prim_data)/sizeof(double);
  return memory;
}
}
