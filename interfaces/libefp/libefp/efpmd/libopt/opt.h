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

#ifndef LIBOPT_OPT_H
#define LIBOPT_OPT_H

#include <stddef.h>

enum opt_result {
	OPT_RESULT_SUCCESS = 0,
	OPT_RESULT_ERROR
};

struct opt_state;

typedef double (*opt_func_t)(size_t, const double *, double *, void *);

struct opt_state *opt_create(size_t);
enum opt_result opt_init(struct opt_state *, size_t, const double *);
void opt_set_func(struct opt_state *, opt_func_t);
void opt_set_user_data(struct opt_state *, void *);
void opt_set_bound(struct opt_state *, size_t, const int *, const double *, const double *);
enum opt_result opt_step(struct opt_state *);
double opt_get_fx(struct opt_state *);
void opt_get_x(struct opt_state *, size_t, double *);
void opt_get_gx(struct opt_state *, size_t, double *);
void opt_shutdown(struct opt_state *);

#endif /* LIBOPT_OPT_H */
