/*-
 * Copyright (c) 2012-2013 Ilya Kaliman
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 */

#ifndef LIBEFP_INT_SHIFT_H
#define LIBEFP_INT_SHIFT_H

static const int shift_ss_x[] = { 0 };
static const int shift_ss_y[] = { 0 };
static const int shift_ss_z[] = { 0 };

static const int shift_sl_x[] = { 0, 1, 0, 0 };
static const int shift_sl_y[] = { 0, 0, 1, 0 };
static const int shift_sl_z[] = { 0, 0, 0, 1 };

static const int shift_sp_x[] = { 1, 0, 0 };
static const int shift_sp_y[] = { 0, 1, 0 };
static const int shift_sp_z[] = { 0, 0, 1 };

static const int shift_sd_x[] = { 2, 0, 0, 1, 1, 0 };
static const int shift_sd_y[] = { 0, 2, 0, 1, 0, 1 };
static const int shift_sd_z[] = { 0, 0, 2, 0, 1, 1 };

static const int shift_sf_x[] = { 3, 0, 0, 2, 2, 1, 0, 1, 0, 1 };
static const int shift_sf_y[] = { 0, 3, 0, 1, 0, 2, 2, 0, 1, 1 };
static const int shift_sf_z[] = { 0, 0, 3, 0, 1, 0, 1, 2, 2, 1 };

static const int shift_ls_x[] = { 0, 5, 0, 0 };
static const int shift_ls_y[] = { 0, 0, 5, 0 };
static const int shift_ls_z[] = { 0, 0, 0, 5 };

static const int shift_ll_x[] = { 0, 1, 0, 0, 5, 6, 5, 5, 0, 1, 0, 0, 0, 1, 0, 0 };
static const int shift_ll_y[] = { 0, 0, 1, 0, 0, 0, 1, 0, 5, 5, 6, 5, 0, 0, 1, 0 };
static const int shift_ll_z[] = { 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 5, 5, 5, 6 };

static const int shift_lp_x[] = { 1, 0, 0, 6, 5, 5, 1, 0, 0, 1, 0, 0 };
static const int shift_lp_y[] = { 0, 1, 0, 0, 1, 0, 5, 6, 5, 0, 1, 0 };
static const int shift_lp_z[] = { 0, 0, 1, 0, 0, 1, 0, 0, 1, 5, 5, 6 };

static const int shift_ld_x[] = { 2, 0, 0, 1, 1, 0, 7, 5, 5, 6, 6, 5, 2, 0, 0, 1, 1, 0, 2, 0, 0, 1, 1, 0 };
static const int shift_ld_y[] = { 0, 2, 0, 1, 0, 1, 0, 2, 0, 1, 0, 1, 5, 7, 5, 6, 5, 6, 0, 2, 0, 1, 0, 1 };
static const int shift_ld_z[] = { 0, 0, 2, 0, 1, 1, 0, 0, 2, 0, 1, 1, 0, 0, 2, 0, 1, 1, 5, 5, 7, 5, 6, 6 };

static const int shift_lf_x[] = { 3, 0, 0, 2, 2, 1, 0, 1, 0, 1, 8, 5, 5, 7, 7, 6, 5, 6, 5, 6, 3, 0, 0, 2, 2, 1, 0, 1, 0, 1, 3, 0, 0, 2, 2, 1, 0, 1, 0, 1 };
static const int shift_lf_y[] = { 0, 3, 0, 1, 0, 2, 2, 0, 1, 1, 0, 3, 0, 1, 0, 2, 2, 0, 1, 1, 5, 8, 5, 6, 5, 7, 7, 5, 6, 6, 0, 3, 0, 1, 0, 2, 2, 0, 1, 1 };
static const int shift_lf_z[] = { 0, 0, 3, 0, 1, 0, 1, 2, 2, 1, 0, 0, 3, 0, 1, 0, 1, 2, 2, 1, 0, 0, 3, 0, 1, 0, 1, 2, 2, 1, 5, 5, 8, 5, 6, 5, 6, 7, 7, 6 };

static const int shift_ps_x[] = { 5, 0, 0 };
static const int shift_ps_y[] = { 0, 5, 0 };
static const int shift_ps_z[] = { 0, 0, 5 };

static const int shift_pl_x[] = { 5, 6, 5, 5, 0, 1, 0, 0, 0, 1, 0, 0 };
static const int shift_pl_y[] = { 0, 0, 1, 0, 5, 5, 6, 5, 0, 0, 1, 0 };
static const int shift_pl_z[] = { 0, 0, 0, 1, 0, 0, 0, 1, 5, 5, 5, 6 };

static const int shift_pp_x[] = { 6, 5, 5, 1, 0, 0, 1, 0, 0 };
static const int shift_pp_y[] = { 0, 1, 0, 5, 6, 5, 0, 1, 0 };
static const int shift_pp_z[] = { 0, 0, 1, 0, 0, 1, 5, 5, 6 };

static const int shift_pd_x[] = { 7, 5, 5, 6, 6, 5, 2, 0, 0, 1, 1, 0, 2, 0, 0, 1, 1, 0 };
static const int shift_pd_y[] = { 0, 2, 0, 1, 0, 1, 5, 7, 5, 6, 5, 6, 0, 2, 0, 1, 0, 1 };
static const int shift_pd_z[] = { 0, 0, 2, 0, 1, 1, 0, 0, 2, 0, 1, 1, 5, 5, 7, 5, 6, 6 };

static const int shift_pf_x[] = { 8, 5, 5, 7, 7, 6, 5, 6, 5, 6, 3, 0, 0, 2, 2, 1, 0, 1, 0, 1, 3, 0, 0, 2, 2, 1, 0, 1, 0, 1 };
static const int shift_pf_y[] = { 0, 3, 0, 1, 0, 2, 2, 0, 1, 1, 5, 8, 5, 6, 5, 7, 7, 5, 6, 6, 0, 3, 0, 1, 0, 2, 2, 0, 1, 1 };
static const int shift_pf_z[] = { 0, 0, 3, 0, 1, 0, 1, 2, 2, 1, 0, 0, 3, 0, 1, 0, 1, 2, 2, 1, 5, 5, 8, 5, 6, 5, 6, 7, 7, 6 };

static const int shift_ds_x[] = { 10, 0, 0, 5, 5, 0 };
static const int shift_ds_y[] = { 0, 10, 0, 5, 0, 5 };
static const int shift_ds_z[] = { 0, 0, 10, 0, 5, 5 };

static const int shift_dl_x[] = { 10, 11, 10, 10, 0, 1, 0, 0, 0, 1, 0, 0, 5, 6, 5, 5, 5, 6, 5, 5, 0, 1, 0, 0 };
static const int shift_dl_y[] = { 0, 0, 1, 0, 10, 10, 11, 10, 0, 0, 1, 0, 5, 5, 6, 5, 0, 0, 1, 0, 5, 5, 6, 5 };
static const int shift_dl_z[] = { 0, 0, 0, 1, 0, 0, 0, 1, 10, 10, 10, 11, 0, 0, 0, 1, 5, 5, 5, 6, 5, 5, 5, 6 };

static const int shift_dp_x[] = { 11, 10, 10, 1, 0, 0, 1, 0, 0, 6, 5, 5, 6, 5, 5, 1, 0, 0 };
static const int shift_dp_y[] = { 0, 1, 0, 10, 11, 10, 0, 1, 0, 5, 6, 5, 0, 1, 0, 5, 6, 5 };
static const int shift_dp_z[] = { 0, 0, 1, 0, 0, 1, 10, 10, 11, 0, 0, 1, 5, 5, 6, 5, 5, 6 };

static const int shift_dd_x[] = { 12, 10, 10, 11, 11, 10, 2, 0, 0, 1, 1, 0, 2, 0, 0, 1, 1, 0, 7, 5, 5, 6, 6, 5, 7, 5, 5, 6, 6, 5, 2, 0, 0, 1, 1, 0 };
static const int shift_dd_y[] = { 0, 2, 0, 1, 0, 1, 10, 12, 10, 11, 10, 11, 0, 2, 0, 1, 0, 1, 5, 7, 5, 6, 5, 6, 0, 2, 0, 1, 0, 1, 5, 7, 5, 6, 5, 6 };
static const int shift_dd_z[] = { 0, 0, 2, 0, 1, 1, 0, 0, 2, 0, 1, 1, 10, 10, 12, 10, 11, 11, 0, 0, 2, 0, 1, 1, 5, 5, 7, 5, 6, 6, 5, 5, 7, 5, 6, 6 };

static const int shift_df_x[] = { 13, 10, 10, 12, 12, 11, 10, 11, 10, 11, 3, 0, 0, 2, 2, 1, 0, 1, 0, 1, 3, 0, 0, 2, 2, 1, 0, 1, 0, 1, 8, 5, 5, 7, 7, 6, 5, 6, 5, 6, 8, 5, 5, 7, 7, 6, 5, 6, 5, 6, 3, 0, 0, 2, 2, 1, 0, 1, 0, 1 };
static const int shift_df_y[] = { 0, 3, 0, 1, 0, 2, 2, 0, 1, 1, 10, 13, 10, 11, 10, 12, 12, 10, 11, 11, 0, 3, 0, 1, 0, 2, 2, 0, 1, 1, 5, 8, 5, 6, 5, 7, 7, 5, 6, 6, 0, 3, 0, 1, 0, 2, 2, 0, 1, 1, 5, 8, 5, 6, 5, 7, 7, 5, 6, 6 };
static const int shift_df_z[] = { 0, 0, 3, 0, 1, 0, 1, 2, 2, 1, 0, 0, 3, 0, 1, 0, 1, 2, 2, 1, 10, 10, 13, 10, 11, 10, 11, 12, 12, 11, 0, 0, 3, 0, 1, 0, 1, 2, 2, 1, 5, 5, 8, 5, 6, 5, 6, 7, 7, 6, 5, 5, 8, 5, 6, 5, 6, 7, 7, 6 };

static const int shift_fs_x[] = { 5, 0, 0, 10, 10, 5, 0, 5, 0, 5 };
static const int shift_fs_y[] = { 0, 15, 0, 5, 0, 10, 10, 0, 5, 5 };
static const int shift_fs_z[] = { 0, 0, 15, 0, 5, 0, 5, 10, 10, 5 };

static const int shift_fl_x[] = { 5, 6, 5, 5, 0, 1, 0, 0, 0, 1, 0, 0, 10, 11, 10, 10, 10, 11, 10, 10, 5, 6, 5, 5, 0, 1, 0, 0, 5, 6, 5, 5, 0, 1, 0, 0, 5, 6, 5, 5 };
static const int shift_fl_y[] = { 0, 0, 1, 0, 15, 15, 16, 15, 0, 0, 1, 0, 5, 5, 6, 5, 0, 0, 1, 0, 10, 10, 11, 10, 10, 10, 11, 10, 0, 0, 1, 0, 5, 5, 6, 5, 5, 5, 6, 5 };
static const int shift_fl_z[] = { 0, 0, 0, 1, 0, 0, 0, 1, 15, 15, 15, 16, 0, 0, 0, 1, 5, 5, 5, 6, 0, 0, 0, 1, 5, 5, 5, 6, 10, 10, 10, 11, 10, 10, 10, 11, 5, 5, 5, 6 };

static const int shift_fp_x[] = { 6, 5, 5, 1, 0, 0, 1, 0, 0, 11, 10, 10, 11, 10, 10, 6, 5, 5, 1, 0, 0, 6, 5, 5, 1, 0, 0, 6, 5, 5 };
static const int shift_fp_y[] = { 0, 1, 0, 15, 16, 15, 0, 1, 0, 5, 6, 5, 0, 1, 0, 10, 11, 10, 10, 11, 10, 0, 1, 0, 5, 6, 5, 5, 6, 5 };
static const int shift_fp_z[] = { 0, 0, 1, 0, 0, 1, 15, 15, 16, 0, 0, 1, 5, 5, 6, 0, 0, 1, 5, 5, 6, 10, 10, 11, 10, 10, 11, 5, 5, 6 };

static const int shift_fd_x[] = { 7, 5, 5, 6, 6, 5, 2, 0, 0, 1, 1, 0, 2, 0, 0, 1, 1, 0, 12, 10, 10, 11, 11, 10, 12, 10, 10, 11, 11, 10, 7, 5, 5, 6, 6, 5, 2, 0, 0, 1, 1, 0, 7, 5, 5, 6, 6, 5, 2, 0, 0, 1, 1, 0, 7, 5, 5, 6, 6, 5 };
static const int shift_fd_y[] = { 0, 2, 0, 1, 0, 1, 15, 17, 15, 16, 15, 16, 0, 2, 0, 1, 0, 1, 5, 7, 5, 6, 5, 6, 0, 2, 0, 1, 0, 1, 10, 12, 10, 11, 10, 11, 10, 12, 10, 11, 10, 11, 0, 2, 0, 1, 0, 1, 5, 7, 5, 6, 5, 6, 5, 7, 5, 6, 5, 6 };
static const int shift_fd_z[] = { 0, 0, 2, 0, 1, 1, 0, 0, 2, 0, 1, 1, 15, 15, 17, 15, 16, 16, 0, 0, 2, 0, 1, 1, 5, 5, 7, 5, 6, 6, 0, 0, 2, 0, 1, 1, 5, 5, 7, 5, 6, 6, 10, 10, 12, 10, 11, 11, 10, 10, 12, 10, 11, 11, 5, 5, 7, 5, 6, 6 };

static const int shift_ff_x[] = { 8, 5, 5, 7, 7, 6, 5, 6, 5, 6, 3, 0, 0, 2, 2, 1, 0, 1, 0, 1, 3, 0, 0, 2, 2, 1, 0, 1, 0, 1, 13, 10, 10, 12, 12, 11, 10, 11, 10, 11, 13, 10, 10, 12, 12, 11, 10, 11, 10, 11, 8, 5, 5, 7, 7, 6, 5, 6, 5, 6, 3, 0, 0, 2, 2, 1, 0, 1, 0, 1, 8, 5, 5, 7, 7, 6, 5, 6, 5, 6, 3, 0, 0, 2, 2, 1, 0, 1, 0, 1, 8, 5, 5, 7, 7, 6, 5, 6, 5, 6 };
static const int shift_ff_y[] = { 0, 3, 0, 1, 0, 2, 2, 0, 1, 1, 15, 18, 15, 16, 15, 17, 17, 15, 16, 16, 0, 3, 0, 1, 0, 2, 2, 0, 1, 1, 5, 8, 5, 6, 5, 7, 7, 5, 6, 6, 0, 3, 0, 1, 0, 2, 2, 0, 1, 1, 10, 13, 10, 11, 10, 12, 12, 10, 11, 11, 10, 13, 10, 11, 10, 12, 12, 10, 11, 11, 0, 3, 0, 1, 0, 2, 2, 0, 1, 1, 5, 8, 5, 6, 5, 7, 7, 5, 6, 6, 5, 8, 5, 6, 5, 7, 7, 5, 6, 6 };
static const int shift_ff_z[] = { 0, 0, 3, 0, 1, 0, 1, 2, 2, 1, 0, 0, 3, 0, 1, 0, 1, 2, 2, 1, 15, 15, 18, 15, 16, 15, 16, 17, 17, 16, 0, 0, 3, 0, 1, 0, 1, 2, 2, 1, 5, 5, 8, 5, 6, 5, 6, 7, 7, 6, 0, 0, 3, 0, 1, 0, 1, 2, 2, 1, 5, 5, 8, 5, 6, 5, 6, 7, 7, 6, 10, 10, 13, 10, 11, 10, 11, 12, 12, 11, 10, 10, 13, 10, 11, 10, 11, 12, 12, 11, 5, 5, 8, 5, 6, 5, 6, 7, 7, 6 };

static const int *shift_table_x[] = { shift_ss_x,
				      shift_sl_x,
				      shift_sp_x,
				      shift_sd_x,
				      shift_sf_x,
				      shift_ls_x,
				      shift_ll_x,
				      shift_lp_x,
				      shift_ld_x,
				      shift_lf_x,
				      shift_ps_x,
				      shift_pl_x,
				      shift_pp_x,
				      shift_pd_x,
				      shift_pf_x,
				      shift_ds_x,
				      shift_dl_x,
				      shift_dp_x,
				      shift_dd_x,
				      shift_df_x,
				      shift_fs_x,
				      shift_fl_x,
				      shift_fp_x,
				      shift_fd_x,
				      shift_ff_x };

static const int *shift_table_y[] = { shift_ss_y,
				      shift_sl_y,
				      shift_sp_y,
				      shift_sd_y,
				      shift_sf_y,
				      shift_ls_y,
				      shift_ll_y,
				      shift_lp_y,
				      shift_ld_y,
				      shift_lf_y,
				      shift_ps_y,
				      shift_pl_y,
				      shift_pp_y,
				      shift_pd_y,
				      shift_pf_y,
				      shift_ds_y,
				      shift_dl_y,
				      shift_dp_y,
				      shift_dd_y,
				      shift_df_y,
				      shift_fs_y,
				      shift_fl_y,
				      shift_fp_y,
				      shift_fd_y,
				      shift_ff_y };

static const int *shift_table_z[] = { shift_ss_z,
				      shift_sl_z,
				      shift_sp_z,
				      shift_sd_z,
				      shift_sf_z,
				      shift_ls_z,
				      shift_ll_z,
				      shift_lp_z,
				      shift_ld_z,
				      shift_lf_z,
				      shift_ps_z,
				      shift_pl_z,
				      shift_pp_z,
				      shift_pd_z,
				      shift_pf_z,
				      shift_ds_z,
				      shift_dl_z,
				      shift_dp_z,
				      shift_dd_z,
				      shift_df_z,
				      shift_fs_z,
				      shift_fl_z,
				      shift_fp_z,
				      shift_fd_z,
				      shift_ff_z };

#endif /* LIBEFP_INT_SHIFT_H */
