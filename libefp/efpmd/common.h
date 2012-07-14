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

#ifndef EFPMD_COMMON_H
#define EFPMD_COMMON_H

#include <stdio.h>
#include <string.h>

#include <efp.h>
#include <phys_const.h>

#define streq(a, b) (strcmp(a, b) == 0)
#define strneq(a, b, n) (strncmp(a, b, n) == 0)

#define ARRAY_SIZE(x) (sizeof(x)/sizeof(x[0]))

#define ANGSTROM_TO_BOHR(x) ((x) / BOHR_RADIUS)
#define BOHR_TO_ANGSTROM(x) ((x) * BOHR_RADIUS)

enum coord_type {
	COORD_TYPE_POINTS = 0,
	COORD_TYPE_XYZABC
};

struct config {
	char *run_type;
	char *fraglib_path;
	char *userlib_path;
	double units_factor;
	enum coord_type coord_type;
	struct efp_opts efp_opts;
	char **potential_file_list;
};

struct sys {
	int n_frag;
	char **frag_name;
	double *frag_coord;
	double *gradient;
};

enum efp_result set_coord(struct efp *,
			  const struct config *,
			  struct sys *);

enum efp_result print_geometry(struct efp *);

enum efp_result print_energy(struct efp *);

enum efp_result print_gradient(struct efp *);

int error(const char *, ...);

int lib_error(enum efp_result);

#endif /* EFPMD_COMMON_H */
