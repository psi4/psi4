/*! \file
    \ingroup CCSORT
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libiwl/iwl.h>
#include <libdpd/dpd.h>

namespace psi { namespace ccsort {

void idx_error(const char *message, int p, int q, int r, int s, int pq, int rs,
	       int pq_sym, int rs_sym, FILE *outfile);
void idx_permute_multipass(dpdfile4 *File, int this_bucket, int **bucket_map, 
			   unsigned long int **bucket_offset, 
                           int p, int q, int r, int s, 
			   int perm_pr, int perm_qs, int perm_prqs,
			   double value, FILE *outfile)
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

  /* Go through the allowed permutations --- NB these are Dirac permutations */

  if(bucket_map[p][q] == this_bucket) {

    /* Get the left and right symmetry blocks */
    pq_sym = p_sym^q_sym;
    rs_sym = r_sym^s_sym;

    /* Get the row and column indices and assign the value */
    pq = Params->rowidx[p][q];
    rs = Params->colidx[r][s];
    if((pq >= Params->rowtot[pq_sym]) || (rs >= Params->coltot[rs_sym]))
      idx_error("MP Params_make: pq, rs", p,q,r,s,pq,rs,pq_sym,rs_sym,outfile);

    offset = bucket_offset[this_bucket][pq_sym];
    File->matrix[pq_sym][pq-offset][rs] = value;
  }

  if(perm_pr) {
    if(bucket_map[r][q] == this_bucket) {
      rq_sym = r_sym^q_sym;
      ps_sym = p_sym^s_sym;
      rq = Params->rowidx[r][q];
      ps = Params->colidx[p][s];
      if((rq >= Params->rowtot[rq_sym]) || (ps >= Params->coltot[ps_sym]))
	idx_error("MP Params_make: rq, ps", p,q,r,s,rq,ps,rq_sym,ps_sym,outfile);

      offset = bucket_offset[this_bucket][rq_sym];
      File->matrix[rq_sym][rq-offset][ps] = value;
    }
  }

  if(perm_qs) {
    if(bucket_map[p][s] == this_bucket) {
      ps_sym = p_sym^s_sym;
      rq_sym = r_sym^q_sym;
      ps = Params->rowidx[p][s];
      rq = Params->colidx[r][q];
      if((ps >= Params->rowtot[ps_sym]) || (rq >= Params->coltot[rq_sym]))
	idx_error("MP Params_make: ps, rq", p,q,r,s,ps,rq,ps_sym,rq_sym,outfile);

      offset = bucket_offset[this_bucket][ps_sym];
      File->matrix[ps_sym][ps-offset][rq] = value;
    }
  }

  if(perm_pr && perm_qs) {
    if(bucket_map[r][s] == this_bucket) {
      rs_sym = r_sym^s_sym;
      pq_sym = p_sym^q_sym;
      rs = Params->rowidx[r][s];
      pq = Params->colidx[p][q];
      if((rs >= Params->rowtot[rs_sym]) || (pq >= Params->coltot[pq_sym]))
	idx_error("MP Params_make: rs, pq", p,q,r,s,rs,pq,rs_sym,pq_sym,outfile);

      offset = bucket_offset[this_bucket][rs_sym];
      File->matrix[rs_sym][rs-offset][pq] = value;
    }
  }

  if(perm_prqs) {
    if(bucket_map[q][p] == this_bucket) {
      qp_sym = q_sym^p_sym;
      sr_sym = s_sym^r_sym;
      qp = Params->rowidx[q][p];
      sr = Params->colidx[s][r];
      if((qp >= Params->rowtot[qp_sym]) || (sr >= Params->coltot[sr_sym]))
	idx_error("MP Params_make: qp, sr", p,q,r,s,qp,sr,qp_sym,sr_sym,outfile);

      offset = bucket_offset[this_bucket][qp_sym];
      File->matrix[qp_sym][qp-offset][sr] = value;
    }

    if(perm_pr) {
      if(bucket_map[q][r] == this_bucket) {
	qr_sym = q_sym^r_sym;
	sp_sym = s_sym^p_sym;
	qr = Params->rowidx[q][r];
	sp = Params->colidx[s][p];
	if((qr >= Params->rowtot[qr_sym])||(sp >= Params->coltot[sp_sym]))
	  idx_error("MP Params_make: qr, sp", p,q,r,s,qr,sp,qr_sym,sp_sym,
		    outfile);

	offset = bucket_offset[this_bucket][qr_sym];
	File->matrix[qr_sym][qr-offset][sp] = value;
      }
    }

    if(perm_qs) {
      if(bucket_map[s][p] == this_bucket) {
	sp_sym = s_sym^p_sym;
	qr_sym = q_sym^r_sym;
	sp = Params->rowidx[s][p];
	qr = Params->colidx[q][r];
	if((sp >= Params->rowtot[sp_sym])||(qr >= Params->coltot[qr_sym]))
	  idx_error("MP Params_make: sp, qr", p,q,r,s,sp,qr,sp_sym,qr_sym,
		    outfile);

	offset = bucket_offset[this_bucket][sp_sym];
	File->matrix[sp_sym][sp-offset][qr] = value;
      }
    }
      
    if(perm_pr && perm_qs) {
      if(bucket_map[s][r] == this_bucket) {
	sr_sym = s_sym^r_sym;
	qp_sym = q_sym^p_sym;
	sr = Params->rowidx[s][r];
	qp = Params->colidx[q][p];
	if((sr >= Params->rowtot[sr_sym])||(qp >= Params->coltot[qp_sym]))
	  idx_error("MP Params_make: sr, qp", p,q,r,s,sr,qp,sr_sym,qp_sym,
		    outfile);

	offset = bucket_offset[this_bucket][sr_sym];
	File->matrix[sr_sym][sr-offset][qp] = value;
      }
    }
  }
}

}} // namespace psi::ccsort
