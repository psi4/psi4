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

#include "common.h"

/* current coordinates from efp struct are used */
void compute_energy(struct state *state, bool do_grad)
{
	struct efp_atom *atoms;
	struct efp_energy efp_energy;
	double xyz[3], xyzabc[6], *grad;
	size_t ifrag, nfrag, iatom, natom;
	int itotal;

	check_fail(efp_compute(state->efp, do_grad));
	check_fail(efp_get_energy(state->efp, &efp_energy));
	check_fail(efp_get_frag_count(state->efp, &nfrag));

	if (do_grad) {
		check_fail(efp_get_gradient(state->efp, state->grad));
		check_fail(efp_get_point_charge_gradient(state->efp, state->grad + 6 * nfrag));
	}

	state->energy = efp_energy.total;

	if (state->ff == NULL)
		return;

	for (ifrag = 0, itotal = 0; ifrag < nfrag; ifrag++) {
		check_fail(efp_get_frag_atom_count(state->efp, ifrag, &natom));
		atoms = malloc(natom * sizeof(struct efp_atom));
		check_fail(efp_get_frag_atoms(state->efp, ifrag, natom, atoms));

		for (iatom = 0; iatom < natom; iatom++, itotal++)
			ff_set_atom_xyz(state->ff, itotal, &atoms[iatom].x);

		free(atoms);
	}

	ff_compute(state->ff, do_grad);

	if (do_grad) {
		for (ifrag = 0, itotal = 0, grad = state->grad; ifrag < nfrag; ifrag++, grad += 6) {
			check_fail(efp_get_frag_xyzabc(state->efp, ifrag, xyzabc));
			check_fail(efp_get_frag_atom_count(state->efp, ifrag, &natom));
			atoms = malloc(natom * sizeof(struct efp_atom));
			check_fail(efp_get_frag_atoms(state->efp, ifrag, natom, atoms));

			for (iatom = 0; iatom < natom; iatom++, itotal++) {
				ff_get_atom_gradient(state->ff, itotal, xyz);

				grad[0] += xyz[0];
				grad[1] += xyz[1];
				grad[2] += xyz[2];

				grad[3] += (atoms[iatom].y - xyzabc[1]) * xyz[2] -
					   (atoms[iatom].z - xyzabc[2]) * xyz[1];
				grad[4] += (atoms[iatom].z - xyzabc[2]) * xyz[0] -
					   (atoms[iatom].x - xyzabc[0]) * xyz[2];
				grad[5] += (atoms[iatom].x - xyzabc[0]) * xyz[1] -
					   (atoms[iatom].y - xyzabc[1]) * xyz[0];
			}

			free(atoms);
		}
	}

	state->energy += ff_get_energy(state->ff);
}
