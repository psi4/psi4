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

#ifndef LIBEFP_INT_H
#define LIBEFP_INT_H

#include "math_util.h"

struct shell {
	char type;       /* shell type - S,L,P,D,F */
	size_t n_funcs;  /* number of functions */
	double *coef;    /* function coefficients */
};

struct xr_atom {
	double x, y, z;
	double znuc;
	size_t n_shells;
	struct shell *shells;
};

void efp_st_int(size_t n_atoms_i,
		const struct xr_atom *atoms_i,
		size_t n_atoms_j,
		const struct xr_atom *atoms_j,
		size_t stride,
		double *s,
		double *t);

void efp_st_int_deriv(size_t n_atoms_i,
		      const struct xr_atom *atoms_i,
		      size_t n_atoms_j,
		      const struct xr_atom *atoms_j,
		      const vec_t *com_i,
		      size_t size_i,
		      size_t size_j,
		      six_t *ds,
		      six_t *dt);

#endif /* LIBEFP_INT_H */
