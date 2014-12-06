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

void NORETURN die(const char *format, ...)
{
	va_list args;

	va_start(args, format);
	vfmsg(stderr, format, args);
	va_end(args);
	fmsg(stderr, "\n");
	fflush(stderr);

#ifdef WITH_MPI
	MPI_Finalize();
#endif

	exit(128);
}

void NORETURN error(const char *format, ...)
{
	char buf[4096];
	va_list args;

	va_start(args, format);
	vsnprintf(buf, sizeof(buf), format, args);
	va_end(args);

	die("ERROR: %s", buf);
}

void check_fail(enum efp_result res)
{
	if (res)
		die("LIBEFP: %s", efp_result_to_string(res));
}

void *xmalloc(size_t size)
{
	void *mem = malloc(size);

	if (!mem)
		error("no memory");

	return mem;
}

void *xcalloc(size_t n, size_t size)
{
	void *mem = calloc(n, size);

	if (!mem)
		error("no memory");

	return mem;
}

void *xrealloc(void *ptr, size_t size)
{
	void *mem = realloc(ptr, size);

	if (!mem)
		error("no memory");

	return mem;
}

void print_vec(const double *vec)
{
	msg("%16.8E %16.8E %16.8E", vec[0], vec[1], vec[2]);
}

void print_geometry(struct efp *efp)
{
	size_t n_frags;
	check_fail(efp_get_frag_count(efp, &n_frags));

	msg("    GEOMETRY (ANGSTROMS)\n\n");

	for (size_t i = 0; i < n_frags; i++) {
		size_t n_atoms;
		check_fail(efp_get_frag_atom_count(efp, i, &n_atoms));

		struct efp_atom atoms[n_atoms];
		check_fail(efp_get_frag_atoms(efp, i, n_atoms, atoms));

		for (size_t a = 0; a < n_atoms; a++) {
			double x = atoms[a].x * BOHR_RADIUS;
			double y = atoms[a].y * BOHR_RADIUS;
			double z = atoms[a].z * BOHR_RADIUS;

			msg("%-16s %12.6lf %12.6lf %12.6lf\n", atoms[a].label, x, y, z);
		}
	}

	size_t n_charges;
	check_fail(efp_get_point_charge_count(efp, &n_charges));

	if (n_charges > 0) {
		double xyz[3 * n_charges];
		check_fail(efp_get_point_charge_coordinates(efp, xyz));

		for (size_t i = 0; i < n_charges; i++) {
			char label[32];
			double x = xyz[3 * i + 0] * BOHR_RADIUS;
			double y = xyz[3 * i + 1] * BOHR_RADIUS;
			double z = xyz[3 * i + 2] * BOHR_RADIUS;

			snprintf(label, sizeof(label), "Q%04zu", i + 1);
			msg("%-16s %12.6lf %12.6lf %12.6lf\n", label, x, y, z);
		}
	}

	msg("\n\n");
}

void print_energy(struct state *state)
{
	struct efp_energy energy;

	check_fail(efp_get_energy(state->efp, &energy));

	msg("    ENERGY COMPONENTS (ATOMIC UNITS)\n\n");
	msg("%30s %16.10lf\n", "ELECTROSTATIC ENERGY", energy.electrostatic);
	msg("%30s %16.10lf\n", "POLARIZATION ENERGY", energy.polarization);
	msg("%30s %16.10lf\n", "DISPERSION ENERGY", energy.dispersion);
	msg("%30s %16.10lf\n", "EXCHANGE REPULSION ENERGY", energy.exchange_repulsion);
	msg("%30s %16.10lf\n", "POINT CHARGES ENERGY", energy.electrostatic_point_charges);
	msg("%30s %16.10lf\n", "CHARGE PENETRATION ENERGY", energy.charge_penetration);
	msg("\n");

	if (state->ff) {
		msg("%30s %16.10lf\n", "FORCE-FIELD ENERGY", ff_get_energy(state->ff));
		msg("\n");
	}

	msg("%30s %16.10lf\n", "TOTAL ENERGY", state->energy);
	msg("\n\n");
}

void print_gradient(struct state *state)
{
	size_t n_frags;
	check_fail(efp_get_frag_count(state->efp, &n_frags));

	for (size_t i = 0; i < n_frags; i++) {
		char name[64];
		check_fail(efp_get_frag_name(state->efp, i, sizeof(name), name));

		msg("    GRADIENT ON FRAGMENT %zu (%s)\n", i + 1, name);
		msg("\nFORCE  ");
		print_vec(state->grad + 6 * i);
		msg("\nTORQUE ");
		print_vec(state->grad + 6 * i + 3);
		msg("\n\n");
	}

	size_t n_charges;
	check_fail(efp_get_point_charge_count(state->efp, &n_charges));

	if (n_charges > 0) {
		double cgrad[3 * n_charges];
		check_fail(efp_get_point_charge_gradient(state->efp, cgrad));
		msg("    GRADIENT ON POINT CHARGES\n\n");

		for (size_t i = 0; i < n_charges; i++) {
			msg("Q%04zu  ", i + 1);
			print_vec(cgrad + 3 * i);
			msg("\n");
		}
	}

	msg("\n");
}

void print_fragment(const char *name, const double *xyzabc, const double *vel)
{
	msg("fragment %s\n", name);

	for (size_t i = 0; i < 6; i++)
		msg(" %14.6e", xyzabc[i]);

	if (vel) {
		msg("\nvelocity\n");

		for (size_t i = 0; i < 6; i++)
			msg(" %14.6e", vel[i]);
	}

	msg("\n\n");
}

void print_charge(double q, double x, double y, double z)
{
	msg("charge %5.2g %14.6e %14.6e %14.6e\n", q, x, y, z);
}

void print_vector(size_t len, const double *vec)
{
	static const size_t CPS = 4;

	for (size_t i = 0; i < len; i += CPS) {
		size_t left = len - i > CPS ? CPS : len - i;

		msg("%8zu  ", i + 1);

		for (size_t ii = 0; ii < left; ii++)
			msg("%16.8E", vec[i + ii]);

		msg("\n");
	}

	msg("\n");
}

void print_matrix(size_t rows, size_t cols, const double *mat)
{
	static const size_t CPS = 4;

	for (size_t j = 0; j < cols; j += CPS) {
		size_t left = cols - j > CPS ? CPS : cols - j;

		msg("    ");

		for (size_t jj = 0; jj < left; jj++)
			msg("%16zu", j + jj + 1);

		msg("\n\n");

		for (size_t i = 0; i < rows; i++) {
			msg("%8zu  ", i + 1);

			for (size_t jj = 0; jj < left; jj++)
				msg("%16.8E", mat[i * cols + j + jj]);

			msg("\n");
		}

		msg("\n\n");
	}
}

vec_t box_from_str(const char *str)
{
	vec_t box;

	if (sscanf(str, "%lf %lf %lf", &box.x, &box.y, &box.z) < 3)
		error("incorrect box format");

	vec_scale(&box, 1.0 / BOHR_RADIUS);
	return box;
}
