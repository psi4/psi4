/*-
 * Copyright (c) 2012-2013 Ilya Kaliman
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
	vfprintf(stderr, format, args);
	va_end(args);
	fprintf(stderr, "\n");
	fflush(stderr);

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
		die("LIBEFP ERROR: %s", efp_result_to_string(res));
}

void *xmalloc(size_t size)
{
	void *mem = malloc(size);

	if (!mem)
		error("NO MEMORY");

	return mem;
}

void *xcalloc(size_t n, size_t size)
{
	void *mem = calloc(n, size);

	if (!mem)
		error("NO MEMORY");

	return mem;
}

void *xrealloc(void *ptr, size_t size)
{
	void *mem = realloc(ptr, size);

	if (!mem)
		error("NO MEMORY");

	return mem;
}

void print_geometry(struct efp *efp)
{
	int n_frags;
	check_fail(efp_get_frag_count(efp, &n_frags));

	printf("    GEOMETRY (ANGSTROMS)\n\n");

	for (int i = 0; i < n_frags; i++) {
		int n_atoms;
		check_fail(efp_get_frag_atom_count(efp, i, &n_atoms));

		struct efp_atom atoms[n_atoms];
		check_fail(efp_get_frag_atoms(efp, i, n_atoms, atoms));

		for (int a = 0; a < n_atoms; a++) {
			double x = atoms[a].x * BOHR_RADIUS;
			double y = atoms[a].y * BOHR_RADIUS;
			double z = atoms[a].z * BOHR_RADIUS;

			printf("%-16s %12.6lf %12.6lf %12.6lf\n", atoms[a].label, x, y, z);
		}
	}

	printf("\n\n");
}

void print_energy(struct efp *efp)
{
	struct efp_energy energy;
	check_fail(efp_get_energy(efp, &energy));

	printf("    ENERGY COMPONENTS (ATOMIC UNITS)\n\n");
	printf("%30s %16.10lf\n", "ELECTROSTATIC ENERGY", energy.electrostatic);
	printf("%30s %16.10lf\n", "CHARGE PENETRATION ENERGY", energy.charge_penetration);
	printf("%30s %16.10lf\n", "POLARIZATION ENERGY", energy.polarization);
	printf("%30s %16.10lf\n", "DISPERSION ENERGY", energy.dispersion);
	printf("%30s %16.10lf\n", "EXCHANGE REPULSION ENERGY", energy.exchange_repulsion);
	printf("\n");
	printf("%30s %16.10lf\n", "TOTAL ENERGY", energy.total);
	printf("\n\n");
}

void print_gradient(struct efp *efp)
{
	int n_frags;
	check_fail(efp_get_frag_count(efp, &n_frags));

	double grad[6 * n_frags];
	check_fail(efp_get_gradient(efp, n_frags, grad));

	for (int i = 0; i < n_frags; i++) {
		char name[64];
		check_fail(efp_get_frag_name(efp, i, sizeof(name), name));

		printf("    GRADIENT ON FRAGMENT %d (%s)\n", i + 1, name);
		printf("\nFORCE  ");

		for (int j = 0; j < 3; j++)
			printf(" %16.8E", grad[6 * i + j]);

		printf("\nTORQUE ");

		for (int j = 3; j < 6; j++)
			printf(" %16.8E", grad[6 * i + j]);

		printf("\n\n");
	}

	printf("\n");
}

void print_fragment(const char *name, const double *xyzabc, const double *vel)
{
	printf("fragment %s\n", name);

	for (int i = 0; i < 6; i++)
		printf(" %14.6e", xyzabc[i]);

	if (vel) {
		printf("\nvelocity\n");

		for (int i = 0; i < 6; i++)
			printf(" %14.6e", vel[i]);
	}

	printf("\n\n");
}

void print_vector(int len, const double *vec)
{
	static const int CPS = 4;

	for (int i = 0; i < len; i += CPS) {
		int left = len - i > CPS ? CPS : len - i;

		printf("%8d  ", i + 1);

		for (int ii = 0; ii < left; ii++)
			printf("%16.8E", vec[i + ii]);

		printf("\n");
	}

	printf("\n");
}

void print_matrix(int rows, int cols, const double *mat)
{
	static const int CPS = 4;

	for (int j = 0; j < cols; j += CPS) {
		int left = cols - j > CPS ? CPS : cols - j;

		printf("    ");

		for (int jj = 0; jj < left; jj++)
			printf("%16d", j + jj + 1);

		printf("\n\n");

		for (int i = 0; i < rows; i++) {
			printf("%8d  ", i + 1);

			for (int jj = 0; jj < left; jj++)
				printf("%16.8E", mat[i * cols + j + jj]);

			printf("\n");
		}

		printf("\n\n");
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
