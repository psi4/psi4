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

#ifndef EFPMD_COMMON_H
#define EFPMD_COMMON_H

#include <assert.h>
#include <ctype.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <efp.h>
#include <math_util.h>

#include "cfg.h"
#include "phys.h"

#define NORETURN __attribute__((noreturn))
#define ARRAY_SIZE(x) (sizeof(x)/sizeof(x[0]))

enum run_type {
	RUN_TYPE_SP,
	RUN_TYPE_GRAD,
	RUN_TYPE_HESS,
	RUN_TYPE_OPT,
	RUN_TYPE_MD
};

enum ensemble_type {
	ENSEMBLE_TYPE_NVE,
	ENSEMBLE_TYPE_NVT,
	ENSEMBLE_TYPE_NPT
};

struct frag {
	char *name;
	double coord[12];
	double vel[6];
};

struct sys {
	int n_frags;
	struct frag *frags;
};

void NORETURN die(const char *, ...);
void NORETURN error(const char *, ...);

void *xmalloc(size_t);
void *xcalloc(size_t, size_t);
void *xrealloc(void *, size_t);

void print_geometry(struct efp *);
void print_energy(struct efp *);
void print_gradient(struct efp *);
void print_fragment(const char *, const double *, const double *);
void print_vector(int, const double *);
void print_matrix(int, int, const double *);

void check_fail(enum efp_result);
struct sys *parse_input(struct cfg *, const char *);
vec_t box_from_str(const char *);
int efp_strcasecmp(const char *, const char *);
int efp_strncasecmp(const char *, const char *, size_t);

#endif /* EFPMD_COMMON_H */
