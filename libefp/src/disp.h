/*-
 * Copyright (c) 2012 Ilya Kaliman
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

#ifndef LIBEFP_DISP_H
#define LIBEFP_DISP_H

static const double disp_weights[] = {
	0.72086099022968040154E-02, 0.17697067815034886394E-01,
	0.30660908596251749739E-01, 0.48381293256249884995E-01,
	0.74878830420650517080E-01, 0.11806515901361630228E+00,
	0.19535413832209084204E+00, 0.35055692324483221824E+00,
	0.71577113554429568336E+00, 0.18140975997632396972E+01,
	0.69792344511487082324E+01, 0.83248093882965845391E+02
};

static inline int
disp_damp_overlap_idx(struct efp *efp, int frag_i, int frag_j,
		      int pt_i, int pt_j)
{
	int n_disp = efp->disp_damp_overlap_offset[efp->n_frag];

	int offset_i = efp->disp_damp_overlap_offset[frag_i] + pt_i;
	int offset_j = efp->disp_damp_overlap_offset[frag_j] + pt_j;

	return offset_i * n_disp + offset_j;
}

#endif /* LIBEFP_DISP_H */
