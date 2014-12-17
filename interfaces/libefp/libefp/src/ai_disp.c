/*-
 * Copyright (c) 2012-2014 Ilya Kaliman
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

#include "balance.h"
#include "private.h"

static const double quad_fact[12] = {
	0.72086099022968040154e-02, 0.17697067815034886394e-01, 0.30660908596251749739e-01,
	0.48381293256249884995e-01, 0.74878830420650517080e-01, 0.11806515901361630228e+00,
	0.19535413832209084204e+00, 0.35055692324483221824e+00, 0.71577113554429568336e+00,
	1.81409759976323969729e+00, 6.97923445114870823247e+00, 8.32480938829658453917e+01
};

static const double quad_freq[12] = {
	0.77932702233253671082e-05, 0.22821071773724297874e-03, 0.15211319247778075068e-02,
	0.60833919905855461032e-02, 0.19223967039304198946e-01, 0.54392829363594207533e-01,
	0.14891668800412559598e+00, 0.42134903703482291779e+00, 1.33149401066630080343e+00,
	5.32498192172462030801e+00, 3.54935126637048206534e+01, 1.03935828835455831714e+03
};

static double
quadrature(const mat_t *tensor, size_t i, size_t j, double de)
{
	double sum = 0.0;

	for (int k = 0; k < 12; k++)
		sum += mat_get(tensor + k, i, j) * quad_fact[k] / (de * de + quad_freq[k]);

	return (sum * de);
}

static double
get_dip_int(struct efp *efp, size_t i_occ, size_t i_vir, size_t axis)
{
	size_t idx, size;

	size = efp->n_ai_core + efp->n_ai_act + efp->n_ai_vir;
	idx = axis * size * size + i_occ * size + (efp->n_ai_core + efp->n_ai_act) + i_vir;

	return (efp->ai_dipole_integrals[idx]);
}

static double
compute_ai_disp_pt(struct efp *efp, size_t fr_idx, size_t pt_idx)
{
	struct frag *frag;
	struct dynamic_polarizable_pt *pt;
	double sum;

	frag = efp->frags + fr_idx;
	pt = frag->dynamic_polarizable_pts + pt_idx;
	sum = 0.0;

	for (size_t i_vir = 0; i_vir < efp->n_ai_vir; i_vir++) {
		double e_vir = efp->ai_orbital_energies[efp->n_ai_core + efp->n_ai_act + i_vir];

		for (size_t i_occ = 0; i_occ < efp->n_ai_core + efp->n_ai_act; i_occ++) {
			double e_occ = efp->ai_orbital_energies[i_occ];

			for (size_t i = 0; i < 3; i++) {
				double dipint_i = get_dip_int(efp, i_occ, i_vir, i);

				for (size_t j = 0; j < 3; j++) {
					double dipint_j = get_dip_int(efp, i_occ, i_vir, j);

					sum += dipint_i * dipint_j * quadrature(pt->tensor, i, j, e_vir - e_occ);
				}
			}
		}
	}

	return (-sum / PI);
}

static void
compute_ai_disp_range(struct efp *efp, size_t from, size_t to, void *data)
{
	double energy = 0.0;

	(void)data;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) reduction(+:energy)
#endif
	for (size_t i = from; i < to; i++) {
		size_t n_pt = efp->frags[i].n_dynamic_polarizable_pts;

		for (size_t j = 0; j < n_pt; j++)
			energy += compute_ai_disp_pt(efp, i, j);
	}

	efp->energy.ai_dispersion += energy;
}

enum efp_result
efp_compute_ai_disp(struct efp *efp)
{
	if (!(efp->opts.terms & EFP_TERM_AI_DISP))
		return (EFP_RESULT_SUCCESS);

	if (efp->do_gradient) {
		efp_log("gradient for AI/EFP dispersion is not implemented");
		return (EFP_RESULT_FATAL);
	}

	efp_balance_work(efp, compute_ai_disp_range, NULL);
	efp_allreduce(&efp->energy.ai_dispersion, 1);

	return (EFP_RESULT_SUCCESS);
}
