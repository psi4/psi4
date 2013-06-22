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

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "stream.h"
#include "private.h"

static inline struct frag *
get_last_frag(struct efp *efp)
{
	return efp->lib[efp->n_lib - 1];
}

static int
tok(struct stream *stream, const char *id)
{
	efp_stream_skip_space(stream);

	if (efp_stream_eol(stream))
		return 0;

	if (strncmp(efp_stream_get_ptr(stream), id, strlen(id)) == 0)
		return efp_stream_advance(stream, strlen(id));

	return 0;
}

static inline int
tok_stop(struct stream *stream)
{
	if (tok(stream, "STOP")) {
		efp_stream_next_line(stream);
		return 1;
	}
	return 0;
}

static inline int
tok_end(struct stream *stream)
{
	if (tok(stream, "$END")) {
		efp_stream_next_line(stream);
		return 1;
	}
	return 0;
}

static int
skip_label(struct stream *stream)
{
	efp_stream_skip_space(stream);
	efp_stream_skip_nonspace(stream);

	return efp_stream_current_char(stream) != '\0';
}

static int
tok_label(struct stream *stream, size_t size, char *val)
{
	efp_stream_skip_space(stream);

	if (efp_stream_eol(stream))
		return 0;

	const char *start = efp_stream_get_ptr(stream);
	efp_stream_skip_nonspace(stream);
	size_t len = efp_stream_get_ptr(stream) - start;

	if (len >= size)
		return 0;

	strncpy(val, start, len);
	val[len] = '\0';
	return 1;
}

static int
tok_int(struct stream *stream, int *val)
{
	return efp_stream_parse_int(stream, val);
}

static int
tok_double(struct stream *stream, double *val)
{
	return efp_stream_parse_double(stream, val);
}

static enum efp_result
parse_coordinates(struct efp *efp, struct stream *stream)
{
	struct frag *frag = get_last_frag(efp);

	efp_stream_next_line(stream);

	while (!efp_stream_eof(stream)) {
		if (tok_stop(stream)) {
			if (frag->n_atoms < 1)
				return EFP_RESULT_SYNTAX_ERROR;

			return EFP_RESULT_SUCCESS;
		}

		struct efp_atom atom;

		if (!tok_label(stream, sizeof(atom.label), atom.label) ||
		    !tok_double(stream, &atom.x) ||
		    !tok_double(stream, &atom.y) ||
		    !tok_double(stream, &atom.z) ||
		    !tok_double(stream, &atom.mass) ||
		    !tok_double(stream, &atom.znuc))
			return EFP_RESULT_SYNTAX_ERROR;

		if (!eq(atom.mass, 0.0)) {
			frag->n_atoms++;
			frag->atoms = realloc(frag->atoms,
				frag->n_atoms * sizeof(struct efp_atom));
			if (!frag->atoms)
				return EFP_RESULT_NO_MEMORY;
			memcpy(frag->atoms + frag->n_atoms - 1, &atom,
				sizeof(struct efp_atom));
		}

		frag->n_multipole_pts++;
		frag->multipole_pts = realloc(frag->multipole_pts,
			frag->n_multipole_pts * sizeof(struct multipole_pt));
		if (!frag->multipole_pts)
			return EFP_RESULT_NO_MEMORY;

		struct multipole_pt *last_pt =
			frag->multipole_pts + frag->n_multipole_pts - 1;

		memset(last_pt, 0, sizeof(struct multipole_pt));
		last_pt->x = atom.x, last_pt->y = atom.y, last_pt->z = atom.z;

		efp_stream_next_line(stream);
	}
	return EFP_RESULT_SYNTAX_ERROR;
}

static enum efp_result
parse_monopoles(struct efp *efp, struct stream *stream)
{
	struct frag *frag = get_last_frag(efp);

	if (!frag->multipole_pts)
		return EFP_RESULT_SYNTAX_ERROR;

	efp_stream_next_line(stream);

	for (int i = 0; i < frag->n_multipole_pts; i++) {
		if (!skip_label(stream) ||
		    !tok_double(stream, &frag->multipole_pts[i].monopole) ||
		    !tok_double(stream, NULL))
			return EFP_RESULT_SYNTAX_ERROR;
		efp_stream_next_line(stream);
	}

	if (!tok_stop(stream))
		return EFP_RESULT_SYNTAX_ERROR;

	return EFP_RESULT_SUCCESS;
}

static enum efp_result
parse_dipoles(struct efp *efp, struct stream *stream)
{
	struct frag *frag = get_last_frag(efp);

	if (!frag->multipole_pts)
		return EFP_RESULT_SYNTAX_ERROR;

	efp_stream_next_line(stream);

	for (int i = 0; i < frag->n_multipole_pts; i++) {
		if (!skip_label(stream) ||
		    !tok_double(stream, &frag->multipole_pts[i].dipole.x) ||
		    !tok_double(stream, &frag->multipole_pts[i].dipole.y) ||
		    !tok_double(stream, &frag->multipole_pts[i].dipole.z))
			return EFP_RESULT_SYNTAX_ERROR;
		efp_stream_next_line(stream);
	}

	if (!tok_stop(stream))
		return EFP_RESULT_SYNTAX_ERROR;

	return EFP_RESULT_SUCCESS;
}

static enum efp_result
parse_quadrupoles(struct efp *efp, struct stream *stream)
{
	struct frag *frag = get_last_frag(efp);

	if (!frag->multipole_pts)
		return EFP_RESULT_SYNTAX_ERROR;

	efp_stream_next_line(stream);

	for (int i = 0; i < frag->n_multipole_pts; i++) {
		if (!skip_label(stream))
			return EFP_RESULT_SYNTAX_ERROR;

		double *q = frag->multipole_pts[i].quadrupole;

		for (int j = 0; j < 6; j++)
			if (!tok_double(stream, q + j))
				return EFP_RESULT_SYNTAX_ERROR;

		efp_stream_next_line(stream);
	}

	if (!tok_stop(stream))
		return EFP_RESULT_SYNTAX_ERROR;

	return EFP_RESULT_SUCCESS;
}

static enum efp_result
parse_octupoles(struct efp *efp, struct stream *stream)
{
	struct frag *frag = get_last_frag(efp);

	if (!frag->multipole_pts)
		return EFP_RESULT_SYNTAX_ERROR;

	efp_stream_next_line(stream);

	for (int i = 0; i < frag->n_multipole_pts; i++) {
		if (!skip_label(stream))
			return EFP_RESULT_SYNTAX_ERROR;

		double *o = frag->multipole_pts[i].octupole;

		for (int j = 0; j < 10; j++)
			if (!tok_double(stream, o + j))
				return EFP_RESULT_SYNTAX_ERROR;

		efp_stream_next_line(stream);
	}

	if (!tok_stop(stream))
		return EFP_RESULT_SYNTAX_ERROR;

	return EFP_RESULT_SUCCESS;
}

static enum efp_result
parse_polarizable_pts(struct efp *efp, struct stream *stream)
{
	struct frag *frag = get_last_frag(efp);
	efp_stream_next_line(stream);

	while (!efp_stream_eof(stream)) {
		if (tok_stop(stream))
			return EFP_RESULT_SUCCESS;

		frag->n_polarizable_pts++;

		size_t size = sizeof(struct polarizable_pt);
		frag->polarizable_pts = realloc(frag->polarizable_pts,
			frag->n_polarizable_pts * size);
		if (!frag->polarizable_pts)
			return EFP_RESULT_NO_MEMORY;

		struct polarizable_pt *pt =
			frag->polarizable_pts + frag->n_polarizable_pts - 1;

		if (!efp_stream_advance(stream, 4))
			return EFP_RESULT_SYNTAX_ERROR;

		if (!tok_double(stream, &pt->x) ||
		    !tok_double(stream, &pt->y) ||
		    !tok_double(stream, &pt->z))
			return EFP_RESULT_SYNTAX_ERROR;

		efp_stream_next_line(stream);
		double m[9];

		for (int i = 0; i < 9; i++)
			if (!tok_double(stream, m + i))
				return EFP_RESULT_SYNTAX_ERROR;

		pt->tensor.xx = m[0];
		pt->tensor.yy = m[1];
		pt->tensor.zz = m[2];
		pt->tensor.xy = m[3];
		pt->tensor.xz = m[4];
		pt->tensor.yz = m[5];
		pt->tensor.yx = m[6];
		pt->tensor.zx = m[7];
		pt->tensor.zy = m[8];

		efp_stream_next_line(stream);
	}

	return EFP_RESULT_SYNTAX_ERROR;
}

static enum efp_result
parse_dynamic_polarizable_pts(struct efp *efp, struct stream *stream)
{
	struct frag *frag = get_last_frag(efp);
	efp_stream_next_line(stream);

	while (!efp_stream_eof(stream)) {
		frag->n_dynamic_polarizable_pts++;

		size_t size = sizeof(struct dynamic_polarizable_pt);
		frag->dynamic_polarizable_pts = realloc(
			frag->dynamic_polarizable_pts,
			frag->n_dynamic_polarizable_pts * size);
		if (!frag->dynamic_polarizable_pts)
			return EFP_RESULT_NO_MEMORY;

		struct dynamic_polarizable_pt *pt =
			frag->dynamic_polarizable_pts +
			frag->n_dynamic_polarizable_pts - 1;

		if (!efp_stream_advance(stream, 5))
			return EFP_RESULT_SYNTAX_ERROR;

		if (!tok_double(stream, &pt->x) ||
		    !tok_double(stream, &pt->y) ||
		    !tok_double(stream, &pt->z))
			return EFP_RESULT_SYNTAX_ERROR;

		efp_stream_next_line(stream);
		double m[9];

		for (int j = 0; j < 9; j++)
			if (!tok_double(stream, m + j))
				return EFP_RESULT_SYNTAX_ERROR;

		pt->trace[0] = (m[0] + m[1] + m[2]) / 3.0;

		efp_stream_next_line(stream);

		if (efp_stream_eof(stream))
			return EFP_RESULT_SYNTAX_ERROR;

		if (strstr(efp_stream_get_ptr(stream), "FOR"))
			break;
	}

	if (efp_stream_eof(stream))
		return EFP_RESULT_SYNTAX_ERROR;

	for (int w = 1; w < 12; w++) {
		for (int i = 0; i < frag->n_dynamic_polarizable_pts; i++) {
			struct dynamic_polarizable_pt *pt =
				frag->dynamic_polarizable_pts + i;

			if (!efp_stream_advance(stream, 5))
				return EFP_RESULT_SYNTAX_ERROR;

			if (!tok_double(stream, &pt->x) ||
			    !tok_double(stream, &pt->y) ||
			    !tok_double(stream, &pt->z))
				return EFP_RESULT_SYNTAX_ERROR;

			efp_stream_next_line(stream);
			double m[9];

			for (int j = 0; j < 9; j++)
				if (!tok_double(stream, m + j))
					return EFP_RESULT_SYNTAX_ERROR;

			pt->trace[w] = (m[0] + m[1] + m[2]) / 3.0;
			efp_stream_next_line(stream);
		}
	}

	if (!tok_stop(stream))
		return EFP_RESULT_SYNTAX_ERROR;

	return EFP_RESULT_SUCCESS;
}

static enum efp_result
parse_projection_basis(struct efp *efp, struct stream *stream)
{
	struct frag *frag = get_last_frag(efp);
	double x, y, z;

	efp_stream_next_line(stream);

	while (!efp_stream_eof(stream)) {
		if (tok_stop(stream))
			return EFP_RESULT_SUCCESS;

		if (!efp_stream_advance(stream, 8))
			return EFP_RESULT_SYNTAX_ERROR;

		if (!tok_double(stream, &x) ||
		    !tok_double(stream, &y) ||
		    !tok_double(stream, &z))
			return EFP_RESULT_SYNTAX_ERROR;

		efp_stream_next_line(stream);
shell:
		if (efp_stream_eof(stream))
			return EFP_RESULT_SYNTAX_ERROR;

		efp_stream_skip_space(stream);

		if (efp_stream_eol(stream)) {
			efp_stream_next_line(stream);
			continue;
		}

		frag->n_xr_shells++;
		frag->xr_shells = realloc(frag->xr_shells,
				frag->n_xr_shells * sizeof(struct shell));

		struct shell *shell = frag->xr_shells + frag->n_xr_shells - 1;

		shell->x = x, shell->y = y, shell->z = z;
		shell->type = efp_stream_get_char(stream);

		if (!strchr("SLPDF", shell->type))
			return EFP_RESULT_SYNTAX_ERROR;

		if (!tok_int(stream, &shell->n_funcs))
			return EFP_RESULT_SYNTAX_ERROR;

		efp_stream_next_line(stream);

		shell->coef = malloc((shell->type == 'L' ? 3 : 2) *
					shell->n_funcs * sizeof(double));

		double *ptr = shell->coef;

		for (int i = 0; i < shell->n_funcs; i++) {
			if (!tok_int(stream, NULL) ||
			    !tok_double(stream, ptr++) ||
			    !tok_double(stream, ptr++))
				return EFP_RESULT_SYNTAX_ERROR;

			if (shell->type == 'L')
				if (!tok_double(stream, ptr++))
					return EFP_RESULT_SYNTAX_ERROR;

			efp_stream_next_line(stream);
		}
		goto shell;
	}

	return EFP_RESULT_SYNTAX_ERROR;
}

static enum efp_result
parse_multiplicity(struct efp *efp, struct stream *stream)
{
	struct frag *frag = get_last_frag(efp);

	if (!tok_int(stream, &frag->multiplicity))
		return EFP_RESULT_SYNTAX_ERROR;

	efp_stream_next_line(stream);

	if (!tok_stop(stream))
		return EFP_RESULT_SYNTAX_ERROR;

	return EFP_RESULT_SUCCESS;
}

static enum efp_result
parse_projection_wf(struct efp *efp, struct stream *stream)
{
	struct frag *frag = get_last_frag(efp);

	if (!tok_int(stream, &frag->n_lmo) ||
	    !tok_int(stream, &frag->xr_wf_size))
		return EFP_RESULT_SYNTAX_ERROR;

	frag->xr_wf = malloc(frag->n_lmo * frag->xr_wf_size * sizeof(double));
	if (!frag->xr_wf)
		return EFP_RESULT_NO_MEMORY;

	efp_stream_next_line(stream);
	double *ptr = frag->xr_wf;

	for (int j = 0; j < frag->n_lmo; j++) {
		for (int i = 0; i < frag->xr_wf_size / 5; i++) {
			if (!efp_stream_advance(stream, 5))
				return EFP_RESULT_SYNTAX_ERROR;

			for (int k = 0; k < 5; k++)
				if (!tok_double(stream, ptr++))
					return EFP_RESULT_SYNTAX_ERROR;

			efp_stream_next_line(stream);
		}

		if (frag->xr_wf_size % 5 == 0)
			continue;

		if (!efp_stream_advance(stream, 5))
			return EFP_RESULT_SYNTAX_ERROR;

		for (int k = 0; k < frag->xr_wf_size % 5; k++)
			if (!tok_double(stream, ptr++))
				return EFP_RESULT_SYNTAX_ERROR;

		efp_stream_next_line(stream);
	}

	return EFP_RESULT_SUCCESS;
}

static enum efp_result
parse_fock_mat(struct efp *efp, struct stream *stream)
{
	struct frag *frag = get_last_frag(efp);
	efp_stream_next_line(stream);

	int size = frag->n_lmo * (frag->n_lmo + 1) / 2;
	frag->xr_fock_mat = malloc(size * sizeof(double));
	if (!frag->xr_fock_mat)
		return EFP_RESULT_NO_MEMORY;

	for (int i = 0; i < size; i++)
		if (!tok_double(stream, frag->xr_fock_mat + i))
			return EFP_RESULT_SYNTAX_ERROR;

	efp_stream_next_line(stream);
	return EFP_RESULT_SUCCESS;
}

static enum efp_result
parse_lmo_centroids(struct efp *efp, struct stream *stream)
{
	struct frag *frag = get_last_frag(efp);
	efp_stream_next_line(stream);

	frag->lmo_centroids = malloc(frag->n_lmo * sizeof(vec_t));
	if (!frag->lmo_centroids)
		return EFP_RESULT_NO_MEMORY;

	for (int i = 0; i < frag->n_lmo; i++) {
		vec_t *ct = frag->lmo_centroids + i;

		if (!skip_label(stream) ||
		    !tok_double(stream, &ct->x) ||
		    !tok_double(stream, &ct->y) ||
		    !tok_double(stream, &ct->z))
			return EFP_RESULT_SYNTAX_ERROR;

		efp_stream_next_line(stream);
	}

	if (!tok_stop(stream))
		return EFP_RESULT_SYNTAX_ERROR;

	return EFP_RESULT_SUCCESS;
}

static enum efp_result
parse_canonvec(struct efp *efp, struct stream *stream)
{
	(void)efp;

	int wf_size;

	if (!tok_int(stream, NULL) ||
	    !tok_int(stream, &wf_size))
		return EFP_RESULT_SYNTAX_ERROR;

	efp_stream_next_line(stream);

	for (int j = 0; j < wf_size; j++) {
		for (int i = 0; i <= (wf_size - 1) / 5; i++) {
			efp_stream_next_line(stream);
		}
	}

	return EFP_RESULT_SUCCESS;
}

static enum efp_result
parse_canonfok(struct efp *efp, struct stream *stream)
{
	(void)efp;

	efp_stream_next_line(stream);

	while (!efp_stream_eof(stream)) {
		if (tok_stop(stream))
			return EFP_RESULT_SUCCESS;
		efp_stream_next_line(stream);
	}

	return EFP_RESULT_SYNTAX_ERROR;
}

static enum efp_result
parse_screen(struct efp *efp, struct stream *stream)
{
	struct frag *frag = get_last_frag(efp);

	double *scr = malloc(frag->n_multipole_pts * sizeof(double));
	if (!scr)
		return EFP_RESULT_NO_MEMORY;

	char type = efp_stream_get_char(stream);

	if (type == '\0' || isspace(type)) {
		if (frag->ai_screen_params)
			return EFP_RESULT_SYNTAX_ERROR;
		frag->ai_screen_params = scr;
	}
	else if (type == '2') {
		if (frag->screen_params)
			return EFP_RESULT_SYNTAX_ERROR;
		frag->screen_params = scr;
	}
	else {
		return EFP_RESULT_UNSUPPORTED_SCREEN;
	}

	efp_stream_next_line(stream);

	for (int i = 0; i < frag->n_multipole_pts; i++) {
		if (!skip_label(stream) ||
		    !tok_double(stream, NULL) ||
		    !tok_double(stream, scr + i))
			return EFP_RESULT_SYNTAX_ERROR;

		efp_stream_next_line(stream);
	}

	if (!tok_stop(stream))
		return EFP_RESULT_SYNTAX_ERROR;

	return EFP_RESULT_SUCCESS;
}

typedef enum efp_result (*parse_fn)(struct efp *, struct stream *);

static parse_fn
get_parse_fn(struct stream *stream)
{
	static const struct {
		const char *label;
		parse_fn fn;
	} funcs[] = {
		{ "COORDINATES",                parse_coordinates },
		{ "MONOPOLES",                  parse_monopoles },
		{ "DIPOLES",                    parse_dipoles },
		{ "QUADRUPOLES",                parse_quadrupoles },
		{ "OCTUPOLES",                  parse_octupoles },
		{ "POLARIZABLE POINTS",         parse_polarizable_pts },
		{ "DYNAMIC POLARIZABLE POINTS", parse_dynamic_polarizable_pts },
		{ "PROJECTION BASIS SET",       parse_projection_basis },
		{ "MULTIPLICITY",               parse_multiplicity },
		{ "PROJECTION WAVEFUNCTION",    parse_projection_wf },
		{ "FOCK MATRIX ELEMENTS",       parse_fock_mat },
		{ "LMO CENTROIDS",              parse_lmo_centroids },
		{ "CANONVEC",                   parse_canonvec },
		{ "CANONFOK",                   parse_canonfok },
		{ "SCREEN",                     parse_screen }
	};

	for (size_t i = 0; i < ARRAY_SIZE(funcs); i++)
		if (tok(stream, funcs[i].label))
			return funcs[i].fn;

	return NULL;
}

static enum efp_result
parse_fragment(struct efp *efp, struct stream *stream)
{
	enum efp_result res;

	while (!efp_stream_eof(stream)) {
		parse_fn fn = get_parse_fn(stream);

		if (!fn) {
			if (tok_end(stream))
				return EFP_RESULT_SUCCESS;

			return EFP_RESULT_SYNTAX_ERROR;
		}

		if ((res = fn(efp, stream)))
			return res;
	}

	return EFP_RESULT_SYNTAX_ERROR;
}

static enum efp_result
parse_file(struct efp *efp, struct stream *stream)
{
	char name[32];
	enum efp_result res;

	while (!efp_stream_eof(stream)) {
		if (efp_stream_get_char(stream) == '\0' ||
		    efp_stream_get_char(stream) != '$') {
			efp_stream_next_line(stream);
			continue;
		}

		if (!tok_label(stream, sizeof(name), name))
			return EFP_RESULT_SYNTAX_ERROR;

		if (efp_find_lib(efp, name))
			return EFP_RESULT_DUPLICATE_PARAMETERS;

		struct frag *frag = calloc(1, sizeof(struct frag));
		if (!frag)
			return EFP_RESULT_NO_MEMORY;

		efp->n_lib++;
		efp->lib = realloc(efp->lib, efp->n_lib * sizeof(struct frag *));
		if (!efp->lib)
			return EFP_RESULT_NO_MEMORY;

		frag->lib = frag;
		strcpy(frag->name, name);
		efp->lib[efp->n_lib - 1] = frag;

		efp_stream_next_line(stream);
		efp_stream_next_line(stream);

		if ((res = parse_fragment(efp, stream)))
			return res;
	}

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_add_potential(struct efp *efp, const char *path)
{
	if (!efp)
		return EFP_RESULT_NOT_INITIALIZED;

	if (!path)
		return EFP_RESULT_INVALID_ARGUMENT;

	enum efp_result res;
	struct stream *stream;

	if ((stream = efp_stream_open(path)) == NULL)
		return EFP_RESULT_FILE_NOT_FOUND;

	efp_stream_set_split_char(stream, '>');
	efp_stream_next_line(stream);
	res = parse_file(efp, stream);
	efp_stream_close(stream);

	return res;
}
