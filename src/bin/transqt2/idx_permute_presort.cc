/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

/*! \file
    \ingroup TRANSQT2
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libiwl/iwl.h>
#include <libdpd/dpd.h>

namespace psi {
  namespace transqt2 {

void idx_error(const char *message, int p, int q, int r, int s, int pq, int rs,
	       int pq_sym, int rs_sym, std::string OutFileRMR);
void idx_permute_presort(dpdfile4 *File, int this_bucket, int **bucket_map, 
			 int **bucket_offset, int p, int q, int r, int s, 
			 double value, std::string OutFileRMR)
{
  int p_sym, q_sym, r_sym, s_sym;
  int pq_sym, rs_sym, rq_sym, ps_sym, qp_sym, sp_sym, sr_sym, qr_sym;
  int pq, rs, rq, ps, qp, sr, qr, sp;
  int perm_pq, perm_rs;
  dpdparams4 *Params;
  int offset;

  Params = File->params;
  perm_pq = Params->perm_pq;
  perm_rs = Params->perm_rs;
  
  /* Get the orbital symmetries */
  p_sym = Params->psym[p]; q_sym = Params->qsym[q];
  r_sym = Params->rsym[r]; s_sym = Params->ssym[s];
  pq_sym = p_sym^q_sym;
  rs_sym = r_sym^s_sym;

  /* The allowed (Mulliken) permutations are very simple in this case */

  if(bucket_map[p][q] == this_bucket) {

    /* Get the row and column indices and assign the value */
    pq = Params->rowidx[p][q];
    rs = Params->colidx[r][s];
    if((pq >= Params->rowtot[pq_sym]) || (rs >= Params->coltot[rs_sym]))
      idx_error("MP Params_make: pq, rs", p,q,r,s,pq,rs,pq_sym,rs_sym,"outfile");

    //printf("p = %d, q = %d, pq_sym = %d, offset = %d, rs = %d, pq = %d\n", p,q,pq_sym, offset, rs, pq);

    offset = bucket_offset[this_bucket][pq_sym];
    File->matrix[pq_sym][pq-offset][rs] = value;
  }

  if(bucket_map[r][s] == this_bucket) {

    rs = Params->rowidx[r][s];
    pq = Params->colidx[p][q];
    if((rs >= Params->rowtot[rs_sym])||(pq >= Params->coltot[pq_sym]))
      idx_error("MP Params_make: rs, pq", p,q,r,s,rs,pq,rs_sym,pq_sym,
		"outfile");

    offset = bucket_offset[this_bucket][rs_sym];
    File->matrix[rs_sym][rs-offset][pq] = value;
  }
}

  } // namespace transqt2
} // namespace psi
