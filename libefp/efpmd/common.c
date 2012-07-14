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

#include <assert.h>
#include <stdarg.h>

#include "common.h"

int error(const char *format, ...)
{
	va_list args;

	fprintf(stderr, "ERROR: ");
	va_start(args, format);
	vfprintf(stderr, format, args);
	va_end(args);
	fprintf(stderr, "\n");
	fflush(stderr);

	return -1;
}

int lib_error(enum efp_result res)
{
	fprintf(stderr, "LIBEFP ");
	return error("%s.", efp_result_to_string(res));
}

enum efp_result set_coord(struct efp *efp,
			  const struct config *config,
			  struct sys *sys)
{
	switch (config->coord_type) {
		case COORD_TYPE_XYZABC:
			return efp_set_coordinates(efp, sys->frag_coord);
		case COORD_TYPE_POINTS:
			return efp_set_coordinates_2(efp, sys->frag_coord);
	}
	assert(0);
}

enum efp_result print_geometry(struct efp *efp)
{
	int n_frag;
	enum efp_result res;

	if ((res = efp_get_frag_count(efp, &n_frag)))
		return res;

	printf("SYSTEM GEOMETRY (ANGSTROM):\n\n");

	for (int i = 0; i < n_frag; i++) {
		int n_atoms;
		if ((res = efp_get_frag_atom_count(efp, i, &n_atoms)))
			return res;

		struct efp_atom atoms[n_atoms];
		if ((res = efp_get_frag_atoms(efp, i, n_atoms, atoms)))
			return res;

		for (int a = 0; a < n_atoms; a++) {
			struct efp_atom *atom = atoms + a;
			printf("%s %12.8lf %12.8lf %12.8lf\n", atom->label,
					BOHR_TO_ANGSTROM(atom->x),
					BOHR_TO_ANGSTROM(atom->y),
					BOHR_TO_ANGSTROM(atom->z));
		}
	}

	printf("\n\n");
	fflush(stdout);
	return EFP_RESULT_SUCCESS;
}

enum efp_result print_energy(struct efp *efp)
{
	enum efp_result res;
	struct efp_energy energy;

	if ((res = efp_get_energy(efp, &energy)))
		return res;

printf("         ELECTROSTATIC ENERGY = %16.10lf\n", energy.electrostatic);
printf("          POLARIZATION ENERGY = %16.10lf\n", energy.polarization);
printf("            DISPERSION ENERGY = %16.10lf\n", energy.dispersion);
printf("    EXCHANGE REPULSION ENERGY = %16.10lf\n", energy.exchange_repulsion);
printf("    CHARGE PENETRATION ENERGY = %16.10lf\n", energy.charge_penetration);
printf("------------------------------------------------\n");
printf("                 TOTAL ENERGY = %16.10lf\n", energy.total);
printf("\n\n");

	fflush(stdout);
	return EFP_RESULT_SUCCESS;
}

enum efp_result print_gradient(struct efp *efp)
{
	int n_frag;
	enum efp_result res;

	if ((res = efp_get_frag_count(efp, &n_frag)))
		return res;

	char frag_name[64];
	double grad[6 * n_frag];

	if ((res = efp_get_gradient(efp, 6 * n_frag, grad)))
		return res;

	for (int i = 0; i < n_frag; i++) {
		if ((res = efp_get_frag_name(efp, i, 64, frag_name)))
			return res;

		double *ptr = grad + 6 * i;

printf("    GRADIENT ON FRAGMENT %d (%s):\n", i, frag_name);
printf(" FORCE: %16.10lf %16.10lf %16.10lf\n", ptr[0], ptr[1], ptr[2]);
printf("TORQUE: %16.10lf %16.10lf %16.10lf\n", ptr[3], ptr[4], ptr[5]);
printf("\n");
	}

	printf("\n");
	fflush(stdout);
	return EFP_RESULT_SUCCESS;
}
