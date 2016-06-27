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

void sim_efield(struct state *state);

static void
print_field(size_t frag_idx, const struct efp_atom *atom, const double *field)
{
	msg("FIELD FOR ATOM %s ON FRAGMENT %zu\n",
			atom->label,
			frag_idx);
	msg("    COORD %12.8lf %12.8lf %12.8lf\n",
			BOHR_RADIUS * atom->x,
			BOHR_RADIUS * atom->y,
			BOHR_RADIUS * atom->z);
	msg("    FIELD %12.8lf %12.8lf %12.8lf\n\n",
			field[0],
			field[1],
			field[2]);
}

void
sim_efield(struct state *state)
{
	size_t n_frags;

	msg("ELECTRIC FIELD JOB\n\n\n");

	print_geometry(state->efp);
	compute_energy(state, false);
	print_energy(state);
	check_fail(efp_get_frag_count(state->efp, &n_frags));

	msg("COORDINATES ARE IN ANGSTROMS\n");
	msg("ELECTRIC FIELD IS IN ATOMIC UNITS\n\n");

	for (size_t i = 0; i < n_frags; i++) {
		double field[3];
		struct efp_atom *atoms;
		size_t n_atoms;

		check_fail(efp_get_frag_atom_count(state->efp, i, &n_atoms));
		atoms = malloc(n_atoms * sizeof(struct efp_atom));
		check_fail(efp_get_frag_atoms(state->efp, i, n_atoms, atoms));

		for (size_t j = 0; j < n_atoms; j++) {
			check_fail(efp_get_electric_field(state->efp, i, &atoms[j].x, field));
			print_field(i + 1, atoms + j, field);
		}

		free(atoms);
	}

	msg("ELECTRIC FIELD JOB COMPLETED SUCCESSFULLY\n");
}
