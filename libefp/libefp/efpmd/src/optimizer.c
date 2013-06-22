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

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "optimizer.h"

struct opt_state {
	int n;
	int m;
	double *x;
	double *l;
	double *u;
	int *nbd;
	double f;
	double *g;
	double factr;
	double pgtol;
	double *wa;
	int *iwa;
	char task[60];
	int iprint;
	char csave[60];
	int lsave[4];
	int isave[44];
	double dsave[29];
	opt_func_t func;
	void *data;
};

void setulb_(int *, int *, double *, double *, double *, int *,
	     double *, double *, double *, double *, double *,
	     int *, char *, int *, char *, int *, int *, double *);

static void call_routine(struct opt_state *state)
{
	setulb_(&state->n,
		&state->m,
		 state->x,
		 state->l,
		 state->u,
		 state->nbd,
		&state->f,
		 state->g,
		&state->factr,
		&state->pgtol,
		 state->wa,
		 state->iwa,
		 state->task,
		&state->iprint,
		 state->csave,
		 state->lsave,
		 state->isave,
		 state->dsave);
}

struct opt_state *opt_create(int n)
{
	assert(n > 0);

	struct opt_state *state = calloc(1, sizeof(struct opt_state));
	assert(state);

	state->n = n;
	state->m = 8;

	state->factr = 0.0;
	state->pgtol = 0.0;

	state->iprint = -1;

	state->x = calloc(n, sizeof(double));
	assert(state->x);
	state->g = calloc(n, sizeof(double));
	assert(state->g);
	state->l = calloc(n, sizeof(double));
	assert(state->l);
	state->u = calloc(n, sizeof(double));
	assert(state->u);

	int wa_size = 5 * n + (2 * n + 11 * state->m + 8) * state->m;
	state->wa = calloc(wa_size, sizeof(double));
	assert(state->wa);

	state->nbd = calloc(n, sizeof(int));
	assert(state->nbd);
	state->iwa = calloc(3 * n, sizeof(int));
	assert(state->iwa);

	return state;
}

enum opt_result opt_init(struct opt_state *state, int n, const double *x)
{
	assert(state);
	assert(n == state->n);
	assert(x);

	memcpy(state->x, x, n * sizeof(double));
	memset(state->task, ' ', sizeof(state->task));

	state->task[0] = 'S';
	state->task[1] = 'T';
	state->task[2] = 'A';
	state->task[3] = 'R';
	state->task[4] = 'T';

	call_routine(state);

	if (strncmp(state->task, "FG_START", strlen("FG_START")) != 0)
		return OPT_RESULT_ERROR;

	state->f = state->func(state->n, state->x, state->g, state->data);

	if (isnan(state->f))
		return OPT_RESULT_ERROR;

	return OPT_RESULT_SUCCESS;
}

void opt_set_func(struct opt_state *state, opt_func_t func)
{
	assert(state);
	state->func = func;
}

void opt_set_user_data(struct opt_state *state, void *data)
{
	assert(state);
	state->data = data;
}

void opt_set_bound(struct opt_state *state, int n, const int *nbd,
		const double *l, const double *u)
{
	assert(state);
	assert(state->n == n);

	memcpy(state->nbd, nbd, n * sizeof(int));
	memcpy(state->l, l, n * sizeof(double));
	memcpy(state->u, u, n * sizeof(double));
}

enum opt_result opt_step(struct opt_state *state)
{
	assert(state);

next:
	call_routine(state);

	if (strncmp(state->task, "FG", strlen("FG")) == 0) {
		state->f = state->func(state->n, state->x, state->g, state->data);

		if (isnan(state->f))
			return OPT_RESULT_ERROR;

		goto next;
	}

	if (strncmp(state->task, "NEW_X", strlen("NEW_X")) == 0)
		return OPT_RESULT_SUCCESS;

	return OPT_RESULT_ERROR;
}

double opt_get_fx(struct opt_state *state)
{
	assert(state);
	return state->f;
}

void opt_get_x(struct opt_state *state, int size, double *out)
{
	assert(state);
	assert(size >= state->n);
	assert(out);

	memcpy(out, state->x, state->n * sizeof(double));
}

void opt_get_gx(struct opt_state *state, int size, double *out)
{
	assert(state);
	assert(size >= state->n);
	assert(out);

	memcpy(out, state->g, state->n * sizeof(double));
}

void opt_shutdown(struct opt_state *state)
{
	if (state) {
		free(state->x);
		free(state->l);
		free(state->u);
		free(state->nbd);
		free(state->g);
		free(state->wa);
		free(state->iwa);
		free(state);
	}
}
