/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include "psi4/libiwl/iwl.h"
#include "psi4/libdpd/dpd.h"

namespace psi { namespace ccdensity {

void idx_error(const char *message, int p, int q, int r, int s, int pq, int rs,
	       int pq_sym, int rs_sym, std::string OutFileRMR);

void idx_permute(dpdfile4 *File, struct iwlbuf *OutBuf,
		 int **bucket_map, int p, int q, int r, int s,
		 int perm_pr, int perm_qs, int perm_prqs,
		 double value, std::string)
{
  int p_sym, q_sym, r_sym, s_sym;
  int pq_sym, rs_sym, rq_sym, ps_sym, qp_sym, sp_sym, sr_sym, qr_sym;
  int pq, rs, rq, ps, qp, sr, qr, sp;
  int perm_pq, perm_rs;
  dpdparams4 *Params;
  int this_bucket;

  Params = File->params;
  perm_pq = Params->perm_pq;
  perm_rs = Params->perm_rs;

  /* Get the orbital symmetries */
  p_sym = Params->psym[p]; q_sym = Params->qsym[q];
  r_sym = Params->rsym[r]; s_sym = Params->ssym[s];

  /* Go through the allowed permutations --- NB these are Dirac permutations */

  /* Get the left and right symmetry blocks */
  pq_sym = p_sym^q_sym;
  rs_sym = r_sym^s_sym;

  /* Get the row and column indices and assign the value */
  pq = Params->rowidx[p][q];
  rs = Params->colidx[r][s];
  if((pq >= Params->rowtot[pq_sym]) || (rs >= Params->coltot[rs_sym]))
      idx_error("Params_make: pq, rs", p,q,r,s,pq,rs,pq_sym,rs_sym,"outfile");

  this_bucket = bucket_map[p][q];
  iwl_buf_wrt_val(&OutBuf[this_bucket], p, q, r, s, value, 0, "outfile", 0);

  if(perm_pr) {
      rq_sym = r_sym^q_sym;
      ps_sym = p_sym^s_sym;
      rq = Params->rowidx[r][q];
      ps = Params->colidx[p][s];
      if((rq >= Params->rowtot[rq_sym]) || (ps >= Params->coltot[ps_sym]))
	  idx_error("Params_make: rq, ps", p,q,r,s,rq,ps,rq_sym,ps_sym,"outfile");

      this_bucket = bucket_map[r][q];
      iwl_buf_wrt_val(&OutBuf[this_bucket], r, q, p, s, value, 0, "outfile", 0);
    }

  if(perm_qs) {
      ps_sym = p_sym^s_sym;
      rq_sym = r_sym^q_sym;
      ps = Params->rowidx[p][s];
      rq = Params->colidx[r][q];
      if((ps >= Params->rowtot[ps_sym]) || (rq >= Params->coltot[rq_sym]))
	  idx_error("Params_make: ps, rq", p,q,r,s,ps,rq,ps_sym,rq_sym,"outfile");

      this_bucket = bucket_map[p][s];
      iwl_buf_wrt_val(&OutBuf[this_bucket], p, s, r, q, value, 0, "outfile", 0);
    }

  if(perm_pr && perm_qs) {
      rs_sym = r_sym^s_sym;
      pq_sym = p_sym^q_sym;
      rs = Params->rowidx[r][s];
      pq = Params->colidx[p][q];
      if((rs >= Params->rowtot[rs_sym]) || (pq >= Params->coltot[pq_sym]))
	  idx_error("Params_make: rs, pq", p,q,r,s,rs,pq,rs_sym,pq_sym,"outfile");

      this_bucket = bucket_map[r][s];
      iwl_buf_wrt_val(&OutBuf[this_bucket], r, s, p, q, value, 0, "outfile", 0);
    }

  if(perm_prqs) {
      qp_sym = q_sym^p_sym;
      sr_sym = s_sym^r_sym;
      qp = Params->rowidx[q][p];
      sr = Params->colidx[s][r];
      if((qp >= Params->rowtot[qp_sym]) || (sr >= Params->coltot[sr_sym]))
	  idx_error("Params_make: qp, sr", p,q,r,s,qp,sr,qp_sym,sr_sym,"outfile");

      this_bucket = bucket_map[q][p];
      iwl_buf_wrt_val(&OutBuf[this_bucket], q, p, s, r, value, 0, "outfile", 0);


      if(perm_pr) {
	  qr_sym = q_sym^r_sym;
	  sp_sym = s_sym^p_sym;
	  qr = Params->rowidx[q][r];
	  sp = Params->colidx[s][p];
	  if((qr >= Params->rowtot[qr_sym])||(sp >= Params->coltot[sp_sym]))
	      idx_error("Params_make: qr, sp", p,q,r,s,qr,sp,qr_sym,sp_sym,
			       "outfile");

	  this_bucket = bucket_map[q][r];
	  iwl_buf_wrt_val(&OutBuf[this_bucket], q, r, s, p, value, 0, "outfile", 0);
	}

      if(perm_qs) {
	  sp_sym = s_sym^p_sym;
	  qr_sym = q_sym^r_sym;
	  sp = Params->rowidx[s][p];
	  qr = Params->colidx[q][r];
	  if((sp >= Params->rowtot[sp_sym])||(qr >= Params->coltot[qr_sym]))
	      idx_error("Params_make: sp, qr", p,q,r,s,sp,qr,sp_sym,qr_sym,
			       "outfile");

	  this_bucket = bucket_map[s][p];
	  iwl_buf_wrt_val(&OutBuf[this_bucket], s, p, q, r, value, 0, "outfile", 0);
	}

      if(perm_pr && perm_qs) {
	  sr_sym = s_sym^r_sym;
	  qp_sym = q_sym^p_sym;
	  sr = Params->rowidx[s][r];
	  qp = Params->colidx[q][p];
	  if((sr >= Params->rowtot[sr_sym])||(qp >= Params->coltot[qp_sym]))
	      idx_error("Params_make: sr, qp", p,q,r,s,sr,qp,sr_sym,qp_sym,
			       "outfile");

	  this_bucket = bucket_map[s][r];
	  iwl_buf_wrt_val(&OutBuf[this_bucket], s, r, q, p, value, 0, "outfile", 0);
	}
    }
}

}} // namespace psi::ccdensity
