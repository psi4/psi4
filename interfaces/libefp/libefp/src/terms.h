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

#ifndef LIBEFP_TERMS_H
#define LIBEFP_TERMS_H

#include "math_util.h"

struct efp;
struct frag;

double efp_frag_frag_elec(struct efp *, size_t, size_t);
double efp_frag_frag_disp(struct efp *, size_t, size_t, const double *, const six_t *);
void efp_frag_frag_xr(struct efp *, size_t, size_t, double *, six_t *, double *, double *);
enum efp_result efp_compute_pol(struct efp *);
enum efp_result efp_compute_ai_elec(struct efp *);
enum efp_result efp_compute_ai_disp(struct efp *);
enum efp_result efp_compute_pol_energy(struct efp *, double *);
void efp_update_elec(struct frag *);
void efp_update_pol(struct frag *);
void efp_update_disp(struct frag *);
void efp_update_xr(struct frag *);

#endif /* LIBEFP_TERMS_H */
