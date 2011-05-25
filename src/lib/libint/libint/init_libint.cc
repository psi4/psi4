#include <stdio.h>
#include <stdlib.h>
#include "libint.h"
#include "hrr_header.h"

extern "C" {
REALTYPE *(*build_eri[6][6][6][6])(Libint_t *, int);
int libint_stack_size[6];

/* This function initializes a matrix of pointers to routines */
/* for computing ERI classes up to (ns|ns) - the base of the library */

void init_libint_base()
{
  build_eri[0][0][1][0] = hrr_order_00p0;
  build_eri[0][0][2][0] = hrr_order_00d0;
  build_eri[0][0][1][1] = hrr_order_00pp;
  build_eri[0][0][3][0] = hrr_order_00f0;
  build_eri[0][0][2][1] = hrr_order_00dp;
  build_eri[0][0][4][0] = hrr_order_00g0;
  build_eri[0][0][3][1] = hrr_order_00fp;
  build_eri[0][0][2][2] = hrr_order_00dd;
  build_eri[0][0][5][0] = hrr_order_00h0;
  build_eri[0][0][4][1] = hrr_order_00gp;
  build_eri[0][0][3][2] = hrr_order_00fd;
  build_eri[0][0][5][1] = hrr_order_00hp;
  build_eri[0][0][4][2] = hrr_order_00gd;
  build_eri[0][0][3][3] = hrr_order_00ff;
  build_eri[0][0][5][2] = hrr_order_00hd;
  build_eri[0][0][4][3] = hrr_order_00gf;
  build_eri[0][0][5][3] = hrr_order_00hf;
  build_eri[0][0][4][4] = hrr_order_00gg;
  build_eri[0][0][5][4] = hrr_order_00hg;
  build_eri[0][0][5][5] = hrr_order_00hh;
  build_eri[1][0][1][0] = hrr_order_p0p0;
  build_eri[1][0][2][0] = hrr_order_p0d0;
  build_eri[1][0][1][1] = hrr_order_p0pp;
  build_eri[1][0][3][0] = hrr_order_p0f0;
  build_eri[1][0][2][1] = hrr_order_p0dp;
  build_eri[1][0][4][0] = hrr_order_p0g0;
  build_eri[1][0][3][1] = hrr_order_p0fp;
  build_eri[1][0][2][2] = hrr_order_p0dd;
  build_eri[1][0][5][0] = hrr_order_p0h0;
  build_eri[1][0][4][1] = hrr_order_p0gp;
  build_eri[1][0][3][2] = hrr_order_p0fd;
  build_eri[1][0][5][1] = hrr_order_p0hp;
  build_eri[1][0][4][2] = hrr_order_p0gd;
  build_eri[1][0][3][3] = hrr_order_p0ff;
  build_eri[1][0][5][2] = hrr_order_p0hd;
  build_eri[1][0][4][3] = hrr_order_p0gf;
  build_eri[1][0][5][3] = hrr_order_p0hf;
  build_eri[1][0][4][4] = hrr_order_p0gg;
  build_eri[1][0][5][4] = hrr_order_p0hg;
  build_eri[1][0][5][5] = hrr_order_p0hh;
  build_eri[2][0][2][0] = hrr_order_d0d0;
  build_eri[2][0][1][1] = hrr_order_d0pp;
  build_eri[2][0][3][0] = hrr_order_d0f0;
  build_eri[2][0][2][1] = hrr_order_d0dp;
  build_eri[2][0][4][0] = hrr_order_d0g0;
  build_eri[2][0][3][1] = hrr_order_d0fp;
  build_eri[2][0][2][2] = hrr_order_d0dd;
  build_eri[2][0][5][0] = hrr_order_d0h0;
  build_eri[2][0][4][1] = hrr_order_d0gp;
  build_eri[2][0][3][2] = hrr_order_d0fd;
  build_eri[2][0][5][1] = hrr_order_d0hp;
  build_eri[2][0][4][2] = hrr_order_d0gd;
  build_eri[2][0][3][3] = hrr_order_d0ff;
  build_eri[2][0][5][2] = hrr_order_d0hd;
  build_eri[2][0][4][3] = hrr_order_d0gf;
  build_eri[2][0][5][3] = hrr_order_d0hf;
  build_eri[2][0][4][4] = hrr_order_d0gg;
  build_eri[2][0][5][4] = hrr_order_d0hg;
  build_eri[2][0][5][5] = hrr_order_d0hh;
  build_eri[1][1][2][0] = hrr_order_ppd0;
  build_eri[1][1][1][1] = hrr_order_pppp;
  build_eri[1][1][3][0] = hrr_order_ppf0;
  build_eri[1][1][2][1] = hrr_order_ppdp;
  build_eri[1][1][4][0] = hrr_order_ppg0;
  build_eri[1][1][3][1] = hrr_order_ppfp;
  build_eri[1][1][2][2] = hrr_order_ppdd;
  build_eri[1][1][5][0] = hrr_order_pph0;
  build_eri[1][1][4][1] = hrr_order_ppgp;
  build_eri[1][1][3][2] = hrr_order_ppfd;
  build_eri[1][1][5][1] = hrr_order_pphp;
  build_eri[1][1][4][2] = hrr_order_ppgd;
  build_eri[1][1][3][3] = hrr_order_ppff;
  build_eri[1][1][5][2] = hrr_order_pphd;
  build_eri[1][1][4][3] = hrr_order_ppgf;
  build_eri[1][1][5][3] = hrr_order_pphf;
  build_eri[1][1][4][4] = hrr_order_ppgg;
  build_eri[1][1][5][4] = hrr_order_pphg;
  build_eri[1][1][5][5] = hrr_order_pphh;
  build_eri[3][0][3][0] = hrr_order_f0f0;
  build_eri[3][0][2][1] = hrr_order_f0dp;
  build_eri[3][0][4][0] = hrr_order_f0g0;
  build_eri[3][0][3][1] = hrr_order_f0fp;
  build_eri[3][0][2][2] = hrr_order_f0dd;
  build_eri[3][0][5][0] = hrr_order_f0h0;
  build_eri[3][0][4][1] = hrr_order_f0gp;
  build_eri[3][0][3][2] = hrr_order_f0fd;
  build_eri[3][0][5][1] = hrr_order_f0hp;
  build_eri[3][0][4][2] = hrr_order_f0gd;
  build_eri[3][0][3][3] = hrr_order_f0ff;
  build_eri[3][0][5][2] = hrr_order_f0hd;
  build_eri[3][0][4][3] = hrr_order_f0gf;
  build_eri[3][0][5][3] = hrr_order_f0hf;
  build_eri[3][0][4][4] = hrr_order_f0gg;
  build_eri[3][0][5][4] = hrr_order_f0hg;
  build_eri[3][0][5][5] = hrr_order_f0hh;
  build_eri[2][1][3][0] = hrr_order_dpf0;
  build_eri[2][1][2][1] = hrr_order_dpdp;
  build_eri[2][1][4][0] = hrr_order_dpg0;
  build_eri[2][1][3][1] = hrr_order_dpfp;
  build_eri[2][1][2][2] = hrr_order_dpdd;
  build_eri[2][1][5][0] = hrr_order_dph0;
  build_eri[2][1][4][1] = hrr_order_dpgp;
  build_eri[2][1][3][2] = hrr_order_dpfd;
  build_eri[2][1][5][1] = hrr_order_dphp;
  build_eri[2][1][4][2] = hrr_order_dpgd;
  build_eri[2][1][3][3] = hrr_order_dpff;
  build_eri[2][1][5][2] = hrr_order_dphd;
  build_eri[2][1][4][3] = hrr_order_dpgf;
  build_eri[2][1][5][3] = hrr_order_dphf;
  build_eri[2][1][4][4] = hrr_order_dpgg;
  build_eri[2][1][5][4] = hrr_order_dphg;
  build_eri[2][1][5][5] = hrr_order_dphh;
  build_eri[4][0][4][0] = hrr_order_g0g0;
  build_eri[4][0][3][1] = hrr_order_g0fp;
  build_eri[4][0][2][2] = hrr_order_g0dd;
  build_eri[4][0][5][0] = hrr_order_g0h0;
  build_eri[4][0][4][1] = hrr_order_g0gp;
  build_eri[4][0][3][2] = hrr_order_g0fd;
  build_eri[4][0][5][1] = hrr_order_g0hp;
  build_eri[4][0][4][2] = hrr_order_g0gd;
  build_eri[4][0][3][3] = hrr_order_g0ff;
  build_eri[4][0][5][2] = hrr_order_g0hd;
  build_eri[4][0][4][3] = hrr_order_g0gf;
  build_eri[4][0][5][3] = hrr_order_g0hf;
  build_eri[4][0][4][4] = hrr_order_g0gg;
  build_eri[4][0][5][4] = hrr_order_g0hg;
  build_eri[4][0][5][5] = hrr_order_g0hh;
  build_eri[3][1][4][0] = hrr_order_fpg0;
  build_eri[3][1][3][1] = hrr_order_fpfp;
  build_eri[3][1][2][2] = hrr_order_fpdd;
  build_eri[3][1][5][0] = hrr_order_fph0;
  build_eri[3][1][4][1] = hrr_order_fpgp;
  build_eri[3][1][3][2] = hrr_order_fpfd;
  build_eri[3][1][5][1] = hrr_order_fphp;
  build_eri[3][1][4][2] = hrr_order_fpgd;
  build_eri[3][1][3][3] = hrr_order_fpff;
  build_eri[3][1][5][2] = hrr_order_fphd;
  build_eri[3][1][4][3] = hrr_order_fpgf;
  build_eri[3][1][5][3] = hrr_order_fphf;
  build_eri[3][1][4][4] = hrr_order_fpgg;
  build_eri[3][1][5][4] = hrr_order_fphg;
  build_eri[3][1][5][5] = hrr_order_fphh;
  build_eri[2][2][4][0] = hrr_order_ddg0;
  build_eri[2][2][3][1] = hrr_order_ddfp;
  build_eri[2][2][2][2] = hrr_order_dddd;
  build_eri[2][2][5][0] = hrr_order_ddh0;
  build_eri[2][2][4][1] = hrr_order_ddgp;
  build_eri[2][2][3][2] = hrr_order_ddfd;
  build_eri[2][2][5][1] = hrr_order_ddhp;
  build_eri[2][2][4][2] = hrr_order_ddgd;
  build_eri[2][2][3][3] = hrr_order_ddff;
  build_eri[2][2][5][2] = hrr_order_ddhd;
  build_eri[2][2][4][3] = hrr_order_ddgf;
  build_eri[2][2][5][3] = hrr_order_ddhf;
  build_eri[2][2][4][4] = hrr_order_ddgg;
  build_eri[2][2][5][4] = hrr_order_ddhg;
  build_eri[2][2][5][5] = hrr_order_ddhh;
  build_eri[5][0][5][0] = hrr_order_h0h0;
  build_eri[5][0][4][1] = hrr_order_h0gp;
  build_eri[5][0][3][2] = hrr_order_h0fd;
  build_eri[5][0][5][1] = hrr_order_h0hp;
  build_eri[5][0][4][2] = hrr_order_h0gd;
  build_eri[5][0][3][3] = hrr_order_h0ff;
  build_eri[5][0][5][2] = hrr_order_h0hd;
  build_eri[5][0][4][3] = hrr_order_h0gf;
  build_eri[5][0][5][3] = hrr_order_h0hf;
  build_eri[5][0][4][4] = hrr_order_h0gg;
  build_eri[5][0][5][4] = hrr_order_h0hg;
  build_eri[5][0][5][5] = hrr_order_h0hh;
  build_eri[4][1][5][0] = hrr_order_gph0;
  build_eri[4][1][4][1] = hrr_order_gpgp;
  build_eri[4][1][3][2] = hrr_order_gpfd;
  build_eri[4][1][5][1] = hrr_order_gphp;
  build_eri[4][1][4][2] = hrr_order_gpgd;
  build_eri[4][1][3][3] = hrr_order_gpff;
  build_eri[4][1][5][2] = hrr_order_gphd;
  build_eri[4][1][4][3] = hrr_order_gpgf;
  build_eri[4][1][5][3] = hrr_order_gphf;
  build_eri[4][1][4][4] = hrr_order_gpgg;
  build_eri[4][1][5][4] = hrr_order_gphg;
  build_eri[4][1][5][5] = hrr_order_gphh;
  build_eri[3][2][5][0] = hrr_order_fdh0;
  build_eri[3][2][4][1] = hrr_order_fdgp;
  build_eri[3][2][3][2] = hrr_order_fdfd;
  build_eri[3][2][5][1] = hrr_order_fdhp;
  build_eri[3][2][4][2] = hrr_order_fdgd;
  build_eri[3][2][3][3] = hrr_order_fdff;
  build_eri[3][2][5][2] = hrr_order_fdhd;
  build_eri[3][2][4][3] = hrr_order_fdgf;
  build_eri[3][2][5][3] = hrr_order_fdhf;
  build_eri[3][2][4][4] = hrr_order_fdgg;
  build_eri[3][2][5][4] = hrr_order_fdhg;
  build_eri[3][2][5][5] = hrr_order_fdhh;
  build_eri[5][1][5][1] = hrr_order_hphp;
  build_eri[5][1][4][2] = hrr_order_hpgd;
  build_eri[5][1][3][3] = hrr_order_hpff;
  build_eri[5][1][5][2] = hrr_order_hphd;
  build_eri[5][1][4][3] = hrr_order_hpgf;
  build_eri[5][1][5][3] = hrr_order_hphf;
  build_eri[5][1][4][4] = hrr_order_hpgg;
  build_eri[5][1][5][4] = hrr_order_hphg;
  build_eri[5][1][5][5] = hrr_order_hphh;
  build_eri[4][2][5][1] = hrr_order_gdhp;
  build_eri[4][2][4][2] = hrr_order_gdgd;
  build_eri[4][2][3][3] = hrr_order_gdff;
  build_eri[4][2][5][2] = hrr_order_gdhd;
  build_eri[4][2][4][3] = hrr_order_gdgf;
  build_eri[4][2][5][3] = hrr_order_gdhf;
  build_eri[4][2][4][4] = hrr_order_gdgg;
  build_eri[4][2][5][4] = hrr_order_gdhg;
  build_eri[4][2][5][5] = hrr_order_gdhh;
  build_eri[3][3][5][1] = hrr_order_ffhp;
  build_eri[3][3][4][2] = hrr_order_ffgd;
  build_eri[3][3][3][3] = hrr_order_ffff;
  build_eri[3][3][5][2] = hrr_order_ffhd;
  build_eri[3][3][4][3] = hrr_order_ffgf;
  build_eri[3][3][5][3] = hrr_order_ffhf;
  build_eri[3][3][4][4] = hrr_order_ffgg;
  build_eri[3][3][5][4] = hrr_order_ffhg;
  build_eri[3][3][5][5] = hrr_order_ffhh;
  build_eri[5][2][5][2] = hrr_order_hdhd;
  build_eri[5][2][4][3] = hrr_order_hdgf;
  build_eri[5][2][5][3] = hrr_order_hdhf;
  build_eri[5][2][4][4] = hrr_order_hdgg;
  build_eri[5][2][5][4] = hrr_order_hdhg;
  build_eri[5][2][5][5] = hrr_order_hdhh;
  build_eri[4][3][5][2] = hrr_order_gfhd;
  build_eri[4][3][4][3] = hrr_order_gfgf;
  build_eri[4][3][5][3] = hrr_order_gfhf;
  build_eri[4][3][4][4] = hrr_order_gfgg;
  build_eri[4][3][5][4] = hrr_order_gfhg;
  build_eri[4][3][5][5] = hrr_order_gfhh;
  build_eri[5][3][5][3] = hrr_order_hfhf;
  build_eri[5][3][4][4] = hrr_order_hfgg;
  build_eri[5][3][5][4] = hrr_order_hfhg;
  build_eri[5][3][5][5] = hrr_order_hfhh;
  build_eri[4][4][5][3] = hrr_order_gghf;
  build_eri[4][4][4][4] = hrr_order_gggg;
  build_eri[4][4][5][4] = hrr_order_gghg;
  build_eri[4][4][5][5] = hrr_order_gghh;
  build_eri[5][4][5][4] = hrr_order_hghg;
  build_eri[5][4][5][5] = hrr_order_hghh;
  build_eri[5][5][5][5] = hrr_order_hhhh;

  libint_stack_size[0] = 1;
  libint_stack_size[1] = 222;
  libint_stack_size[2] = 3193;
  libint_stack_size[3] = 34956;
  libint_stack_size[4] = 219445;
  libint_stack_size[5] = 759382;
}
/* These functions initialize library objects */
/* Library objects operate independently of each other */
int init_libint(Libint_t *libint, int max_am, int max_num_prim_quartets)
{
  int memory = 0;

  if (max_am >= LIBINT_MAX_AM) return -1;
  libint->int_stack = (REALTYPE *) malloc(libint_stack_size[max_am]*sizeof(REALTYPE));
  memory += libint_stack_size[max_am];
  libint->PrimQuartet = (prim_data *) malloc(max_num_prim_quartets*sizeof(prim_data));
  memory += max_num_prim_quartets*sizeof(prim_data)/sizeof(REALTYPE);
  return memory;
}

void free_libint(Libint_t *libint)
{
  if (libint->int_stack != NULL) {
    free(libint->int_stack);
    libint->int_stack = NULL;
  }
  if (libint->PrimQuartet != NULL) {
    free(libint->PrimQuartet);
    libint->PrimQuartet = NULL;
  }

  return;
}

int libint_storage_required(int max_am, int max_num_prim_quartets)
{
  int memory = 0;

  if (max_am >= LIBINT_MAX_AM) return -1;
  memory += libint_stack_size[max_am];
  memory += max_num_prim_quartets*sizeof(prim_data)/sizeof(REALTYPE);
  return memory;
}

}
