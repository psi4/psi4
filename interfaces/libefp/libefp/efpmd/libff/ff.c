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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ff.h"

#ifdef WITH_TINKER

#define KCALMOL_TO_AU (1.0 / 627.50947)
#define ANGSTROM_TO_BOHR (1.0 / 0.52917721092)

/* must be consistent with tinker */
#define MAXARG 20
#define MAXATM 100000

void initial_(void);
void getxyz_(void);
void mechanic_(void);
double energy_(void);
void gradient_(double *, double *);
void final_(void);

struct atoms {
	double x[MAXATM];
	double y[MAXATM];
	double z[MAXATM];
	int n;
	int type[MAXATM];
};

struct argue {
	int narg;
	int listarg[MAXARG + 1];
	char arg[120 * (MAXARG + 1)];
};

extern struct atoms atoms_;
extern struct argue argue_;

struct ff {
	double energy;
	double grad[3 * MAXATM];
};

struct ff *
ff_create(void)
{
	struct ff *ff;

	ff = calloc(1, sizeof(*ff));

	initial_();

	argue_.narg = 0;
	memset(argue_.arg, ' ', 120 * (MAXARG + 1));

	return (ff);
}

int
ff_load_geometry(struct ff *ff, const char *path)
{
	FILE *fp;
	int i;

	(void)ff;

	if ((fp = fopen(path, "r")) == NULL)
		return (0);

	fclose(fp);

	argue_.narg++;
	argue_.listarg[argue_.narg] = 1;

	for (i = 0; *path; i++, path++)
		argue_.arg[120 * argue_.narg + i] = *path;

	getxyz_();
	return (1);
}

int
ff_load_parameters(struct ff *ff, const char *path)
{
	FILE *fp;
	int i;

	(void)ff;

	if ((fp = fopen(path, "r")) == NULL)
		return (0);

	fclose(fp);

	argue_.narg++;
	argue_.listarg[argue_.narg] = 1;

	for (i = 0; *path; i++, path++)
		argue_.arg[120 * argue_.narg + i] = *path;

	mechanic_();
	return (1);
}

int
ff_get_atom_count(struct ff *ff)
{
	(void)ff;

	return (atoms_.n);
}

void
ff_get_atom_xyz(struct ff *ff, int idx, double *xyz)
{
	(void)ff;

	xyz[0] = atoms_.x[idx] * ANGSTROM_TO_BOHR;
	xyz[1] = atoms_.y[idx] * ANGSTROM_TO_BOHR;
	xyz[2] = atoms_.z[idx] * ANGSTROM_TO_BOHR;
}

void
ff_set_atom_xyz(struct ff *ff, int idx, const double *xyz)
{
	(void)ff;

	atoms_.x[idx] = xyz[0] / ANGSTROM_TO_BOHR;
	atoms_.y[idx] = xyz[1] / ANGSTROM_TO_BOHR;
	atoms_.z[idx] = xyz[2] / ANGSTROM_TO_BOHR;
}

void
ff_compute(struct ff *ff, int do_grad)
{
	int i;

	if (do_grad)
		gradient_(&ff->energy, ff->grad);
	else
		ff->energy = energy_();

	if (do_grad)
		for (i = 0; i < 3 * MAXATM; i++)
			ff->grad[i] *= KCALMOL_TO_AU / ANGSTROM_TO_BOHR;

	ff->energy *= KCALMOL_TO_AU;
}

double
ff_get_energy(struct ff *ff)
{
	return (ff->energy);
}

void
ff_get_atom_gradient(struct ff *ff, int idx, double *grad)
{
	memcpy(grad, ff->grad + 3 * idx, 3 * sizeof(double));
}

void
ff_free(struct ff *ff)
{
	if (ff) {
		free(ff);
		final_();
	}
}

#else /* WITH_TINKER */

struct ff {
	int dummy;
};

struct ff *
ff_create(void)
{
	struct ff *ff;

	ff = calloc(1, sizeof(struct ff));

	return (ff);
}

int
ff_load_geometry(struct ff *ff, const char *path)
{
	(void)ff;
	(void)path;

	return (0);
}

int
ff_load_parameters(struct ff *ff, const char *path)
{
	(void)ff;
	(void)path;

	return (0);
}

int
ff_get_atom_count(struct ff *ff)
{
	(void)ff;

	return (0);
}

void
ff_get_atom_xyz(struct ff *ff, int idx, double *xyz)
{
	(void)ff;
	(void)idx;
	(void)xyz;
}

void
ff_set_atom_xyz(struct ff *ff, int idx, const double *xyz)
{
	(void)ff;
	(void)idx;
	(void)xyz;
}

void
ff_compute(struct ff *ff, int do_grad)
{
	(void)ff;
	(void)do_grad;
}

double
ff_get_energy(struct ff *ff)
{
	(void)ff;

	return (0.0);
}

void
ff_get_atom_gradient(struct ff *ff, int idx, double *xyz)
{
	(void)ff;
	(void)idx;
	(void)xyz;
}

void
ff_free(struct ff *ff)
{
	free(ff);
}

#endif /* WITH_TINKER */
