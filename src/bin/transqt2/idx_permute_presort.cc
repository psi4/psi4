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
	       int pq_sym, int rs_sym, FILE *outfile);
void idx_permute_presort(dpdfile4 *File, int this_bucket, int **bucket_map, 
			 int **bucket_offset, int p, int q, int r, int s, 
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
  pq_sym = p_sym^q_sym;
  rs_sym = r_sym^s_sym;

  /* The allowed (Mulliken) permutations are very simple in this case */

  if(bucket_map[p][q] == this_bucket) {

    /* Get the row and column indices and assign the value */
    pq = Params->rowidx[p][q];
    rs = Params->colidx[r][s];
    if((pq >= Params->rowtot[pq_sym]) || (rs >= Params->coltot[rs_sym]))
      idx_error("MP Params_make: pq, rs", p,q,r,s,pq,rs,pq_sym,rs_sym,outfile);

    //printf("p = %d, q = %d, pq_sym = %d, offset = %d, rs = %d, pq = %d\n", p,q,pq_sym, offset, rs, pq);

    offset = bucket_offset[this_bucket][pq_sym];
    File->matrix[pq_sym][pq-offset][rs] = value;
  }

  if(bucket_map[r][s] == this_bucket) {

    rs = Params->rowidx[r][s];
    pq = Params->colidx[p][q];
    if((rs >= Params->rowtot[rs_sym])||(pq >= Params->coltot[pq_sym]))
      idx_error("MP Params_make: rs, pq", p,q,r,s,rs,pq,rs_sym,pq_sym,
		outfile);

    offset = bucket_offset[this_bucket][rs_sym];
    File->matrix[rs_sym][rs-offset][pq] = value;
  }
}

  } // namespace transqt2
} // namespace psi
