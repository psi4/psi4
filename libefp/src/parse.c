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

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "efp_private.h"
#include "parse.h"

#define streq(a, b) (strcmp((a), (b)) == 0)

struct stream {
	char *buffer;
	char *ptr;
	FILE *in;
};

static inline struct frag *
get_last_frag(struct efp *efp)
{
	return efp->lib + efp->n_lib - 1;
}

static char *
read_line(FILE *in)
{
	int size = 128;
	int i = 0;
	char *buffer = malloc(size);

	if (!buffer)
		return NULL;

	for(;;) {
		int ch = getc(in);

		switch(ch) {
		case EOF:
			if (i == 0) {
				free(buffer);
				return NULL;
			}

		case '\n':
			if (i == size)
				buffer = realloc(buffer, size + 1);
			buffer[i] = '\0';
			return buffer;

		case '>':
			ch = getc(in);

			if (ch == '\n')
				continue;

		default:
			buffer[i++] = ch;

			if (i == size) {
				size *= 2;
				buffer = realloc(buffer, size);
			}
		}
	}
	return NULL;
}

static void
next_line(struct stream *stream)
{
	if (stream->buffer)
		free(stream->buffer);

	stream->buffer = read_line(stream->in);
	stream->ptr = stream->buffer;
}

static void
skip_space(struct stream *stream)
{
	while (*stream->ptr && isspace(*stream->ptr))
		stream->ptr++;
}

static int
tok(struct stream *stream, const char *id)
{
	if (!stream->ptr)
		return 0;

	skip_space(stream);
	size_t len = strlen(id);

	if (strncmp(stream->ptr, id, len) == 0) {
		stream->ptr += len;
		return 1;
	}

	return 0;
}

static inline int
tok_stop(struct stream *stream)
{
	if (tok(stream, "STOP")) {
		next_line(stream);
		return 1;
	}
	return 0;
}

static inline int
tok_end(struct stream *stream)
{
	if (tok(stream, "$END")) {
		next_line(stream);
		return 1;
	}
	return 0;
}

static int
tok_label(struct stream *stream, char **val)
{
	skip_space(stream);
	size_t len = 0;

	while (stream->ptr[len] && !isspace(stream->ptr[len]))
		len++;

	if (val)
		*val = strndup(stream->ptr, len);

	stream->ptr += len;
	return 1;
}

static int
tok_int(struct stream *stream, int *val)
{
	skip_space(stream);

	char *endptr;
	int x = strtol(stream->ptr, &endptr, 10);

	if (endptr == stream->ptr)
		return 0;

	if (val)
		*val = x;

	stream->ptr = endptr;
	return 1;
}

static int
tok_double(struct stream *stream, double *val)
{
	skip_space(stream);

	char *endptr;
	double x = strtod(stream->ptr, &endptr);

	if (endptr == stream->ptr)
		return 0;

	if (val)
		*val = x;

	stream->ptr = endptr;
	return 1;
}

static enum efp_result
parse_coordinates(struct efp *efp, struct stream *stream)
{
	struct frag *frag = get_last_frag(efp);

	next_line(stream);

	while (stream->ptr) {
		if (tok_stop(stream)) {
			if (frag->n_atoms < 1)
				return EFP_RESULT_SYNTAX_ERROR;

			return EFP_RESULT_SUCCESS;
		}

		char *label;
		struct efp_atom atom;

		if (!tok_label(stream, &label) ||
		    !tok_double(stream, &atom.x) ||
		    !tok_double(stream, &atom.y) ||
		    !tok_double(stream, &atom.z) ||
		    !tok_double(stream, &atom.mass) ||
		    !tok_double(stream, &atom.znuc))
			return EFP_RESULT_SYNTAX_ERROR;

		if (strlen(label) >= ARRAY_SIZE(atom.label)) {
			free(label);
			return EFP_RESULT_SYNTAX_ERROR;
		}

		strcpy(atom.label, label);
		free(label);

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

		next_line(stream);
	}
	return EFP_RESULT_SYNTAX_ERROR;
}

static enum efp_result
parse_monopoles(struct efp *efp, struct stream *stream)
{
	struct frag *frag = get_last_frag(efp);

	if (!frag->multipole_pts)
		return EFP_RESULT_SYNTAX_ERROR;

	next_line(stream);

	for (int i = 0; i < frag->n_multipole_pts; i++) {
		if (!tok_label(stream, NULL) ||
		    !tok_double(stream, &frag->multipole_pts[i].monopole) ||
		    !tok_double(stream, NULL))
			return EFP_RESULT_SYNTAX_ERROR;
		next_line(stream);
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

	next_line(stream);

	for (int i = 0; i < frag->n_multipole_pts; i++) {
		if (!tok_label(stream, NULL) ||
		    !tok_double(stream, &frag->multipole_pts[i].dipole.x) ||
		    !tok_double(stream, &frag->multipole_pts[i].dipole.y) ||
		    !tok_double(stream, &frag->multipole_pts[i].dipole.z))
			return EFP_RESULT_SYNTAX_ERROR;
		next_line(stream);
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

	next_line(stream);

	for (int i = 0; i < frag->n_multipole_pts; i++) {
		if (!tok_label(stream, NULL))
			return EFP_RESULT_SYNTAX_ERROR;

		double *q = frag->multipole_pts[i].quadrupole;

		for (int j = 0; j < 6; j++)
			if (!tok_double(stream, q + j))
				return EFP_RESULT_SYNTAX_ERROR;

		next_line(stream);
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

	next_line(stream);

	for (int i = 0; i < frag->n_multipole_pts; i++) {
		if (!tok_label(stream, NULL))
			return EFP_RESULT_SYNTAX_ERROR;

		double *o = frag->multipole_pts[i].octupole;

		for (int j = 0; j < 10; j++)
			if (!tok_double(stream, o + j))
				return EFP_RESULT_SYNTAX_ERROR;

		next_line(stream);
	}

	if (!tok_stop(stream))
		return EFP_RESULT_SYNTAX_ERROR;

	return EFP_RESULT_SUCCESS;
}

static enum efp_result
parse_polarizable_pts(struct efp *efp, struct stream *stream)
{
	struct frag *frag = get_last_frag(efp);
	next_line(stream);

	while (stream->ptr) {
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

		stream->ptr += 4;

		if (!tok_double(stream, &pt->x) ||
		    !tok_double(stream, &pt->y) ||
		    !tok_double(stream, &pt->z))
			return EFP_RESULT_SYNTAX_ERROR;

		next_line(stream);
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

		next_line(stream);
	}
	return EFP_RESULT_SYNTAX_ERROR;
}

static enum efp_result
parse_dynamic_polarizable_pts(struct efp *efp, struct stream *stream)
{
	struct frag *frag = get_last_frag(efp);
	next_line(stream);

	while (stream->ptr) {
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

		stream->ptr += 5;

		if (!tok_double(stream, &pt->x) ||
		    !tok_double(stream, &pt->y) ||
		    !tok_double(stream, &pt->z))
			return EFP_RESULT_SYNTAX_ERROR;

		next_line(stream);
		double m[9];

		for (int j = 0; j < 9; j++)
			if (!tok_double(stream, m + j))
				return EFP_RESULT_SYNTAX_ERROR;

		pt->trace[0] = (m[0] + m[1] + m[2]) / 3.0;

		next_line(stream);
		if (strstr(stream->ptr, "FOR"))
			break;
	}

	if (!stream->ptr)
		return EFP_RESULT_SYNTAX_ERROR;

	for (int w = 1; w < 12; w++) {
		for (int i = 0; i < frag->n_dynamic_polarizable_pts; i++) {
			struct dynamic_polarizable_pt *pt =
				frag->dynamic_polarizable_pts + i;

			stream->ptr += 5;

			if (!tok_double(stream, &pt->x) ||
			    !tok_double(stream, &pt->y) ||
			    !tok_double(stream, &pt->z))
				return EFP_RESULT_SYNTAX_ERROR;

			next_line(stream);
			double m[9];

			for (int j = 0; j < 9; j++)
				if (!tok_double(stream, m + j))
					return EFP_RESULT_SYNTAX_ERROR;

			pt->trace[w] = (m[0] + m[1] + m[2]) / 3.0;
			next_line(stream);
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

	next_line(stream);

	while (stream->ptr) {
		if (tok_stop(stream))
			return EFP_RESULT_SUCCESS;

		stream->ptr += 8;

		if (!tok_double(stream, &x) ||
		    !tok_double(stream, &y) ||
		    !tok_double(stream, &z))
			return EFP_RESULT_SYNTAX_ERROR;

		next_line(stream);
shell:
		if (!stream->ptr)
			return EFP_RESULT_SYNTAX_ERROR;

		skip_space(stream);

		if (!*stream->ptr) {
			next_line(stream);
			continue;
		}

		frag->n_xr_shells++;
		frag->xr_shells = realloc(frag->xr_shells,
				frag->n_xr_shells * sizeof(struct shell));

		struct shell *shell = frag->xr_shells + frag->n_xr_shells - 1;

		shell->x = x, shell->y = y, shell->z = z;
		shell->type = *stream->ptr++;

		if (!strchr("SLPDF", shell->type))
			return EFP_RESULT_SYNTAX_ERROR;

		if (!tok_int(stream, &shell->n_funcs))
			return EFP_RESULT_SYNTAX_ERROR;

		next_line(stream);

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

			next_line(stream);
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

	next_line(stream);

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

	next_line(stream);
	double *ptr = frag->xr_wf;

	for (int j = 0; j < frag->n_lmo; j++) {
		for (int i = 0; i < frag->xr_wf_size / 5; i++) {
			stream->ptr += 5;
			for (int k = 0; k < 5; k++)
				if (!tok_double(stream, ptr++))
					return EFP_RESULT_SYNTAX_ERROR;
			next_line(stream);
		}
		for (int k = 0; k < frag->xr_wf_size % 5; k++) {
			stream->ptr += 5;
			if (!tok_double(stream, ptr++))
				return EFP_RESULT_SYNTAX_ERROR;
		}
		if (frag->xr_wf_size % 5)
			next_line(stream);
	}
	return EFP_RESULT_SUCCESS;
}

static enum efp_result
parse_fock_mat(struct efp *efp, struct stream *stream)
{
	struct frag *frag = get_last_frag(efp);
	next_line(stream);

	int size = frag->n_lmo * (frag->n_lmo + 1) / 2;
	frag->xr_fock_mat = malloc(size * sizeof(double));
	if (!frag->xr_fock_mat)
		return EFP_RESULT_NO_MEMORY;

	for (int i = 0; i < size; i++)
		if (!tok_double(stream, frag->xr_fock_mat + i))
			return EFP_RESULT_SYNTAX_ERROR;

	next_line(stream);
	return EFP_RESULT_SUCCESS;
}

static enum efp_result
parse_lmo_centroids(struct efp *efp, struct stream *stream)
{
	struct frag *frag = get_last_frag(efp);
	next_line(stream);

	frag->lmo_centroids = malloc(frag->n_lmo * sizeof(vec_t));
	if (!frag->lmo_centroids)
		return EFP_RESULT_NO_MEMORY;

	for (int i = 0; i < frag->n_lmo; i++) {
		vec_t *ct = frag->lmo_centroids + i;

		if (!tok_label(stream, NULL) ||
		    !tok_double(stream, &ct->x) ||
		    !tok_double(stream, &ct->y) ||
		    !tok_double(stream, &ct->z))
			return EFP_RESULT_SYNTAX_ERROR;

		next_line(stream);
	}

	if (!tok_stop(stream))
		return EFP_RESULT_SYNTAX_ERROR;

	return EFP_RESULT_SUCCESS;
}

static enum efp_result
parse_canonvec(struct efp *efp, struct stream *stream)
{
	int wf_size;

	if (!tok_int(stream, NULL) ||
	    !tok_int(stream, &wf_size))
		return EFP_RESULT_SYNTAX_ERROR;

	next_line(stream);

	for (int j = 0; j < wf_size; j++) {
		for (int i = 0; i <= (wf_size - 1) / 5; i++) {
			next_line(stream);
		}
	}
	return EFP_RESULT_SUCCESS;
}

static enum efp_result
parse_canonfok(struct efp *efp, struct stream *stream)
{
	next_line(stream);

	while (stream->ptr) {
		if (tok_stop(stream))
			return EFP_RESULT_SUCCESS;
		next_line(stream);
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

	char type = *stream->ptr;

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

	next_line(stream);

	for (int i = 0; i < frag->n_multipole_pts; i++) {
		if (!tok_label(stream, NULL) ||
		    !tok_double(stream, NULL) ||
		    !tok_double(stream, scr + i))
			return EFP_RESULT_SYNTAX_ERROR;

		next_line(stream);
	}

	if (!tok_stop(stream))
		return EFP_RESULT_SYNTAX_ERROR;

	return EFP_RESULT_SUCCESS;
}

typedef enum efp_result (*parse_fn)(struct efp *efp, struct stream *stream);

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

	while (stream->ptr) {
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
	enum efp_result res;

	while (stream->ptr) {
		if (stream->ptr[0] == '\0' || stream->ptr[1] != '$') {
			next_line(stream);
			continue;
		}

		stream->ptr += 2;

		efp->n_lib++;
		efp->lib = realloc(efp->lib, efp->n_lib * sizeof(struct frag));
		if (!efp->lib)
			return EFP_RESULT_NO_MEMORY;

		struct frag *frag = get_last_frag(efp);
		memset(frag, 0, sizeof(struct frag));

		tok_label(stream, &frag->name);
		next_line(stream);
		next_line(stream);

		if ((res = parse_fragment(efp, stream)))
			return res;
	}
	return EFP_RESULT_SUCCESS;
}

static void
set_frag_com(struct frag *frag)
{
	frag->x = 0.0;
	frag->y = 0.0;
	frag->z = 0.0;

	double mass = 0.0;

	for (int i = 0; i < frag->n_atoms; i++) {
		const struct efp_atom *atom = frag->atoms + i;

		frag->x += atom->x * atom->mass;
		frag->y += atom->y * atom->mass;
		frag->z += atom->z * atom->mass;

		mass += atom->mass;
	}

	frag->x /= mass;
	frag->y /= mass;
	frag->z /= mass;
}

enum efp_result
efp_read_potential(struct efp *efp, const char **files)
{
	for (int i = 0; files[i]; i++) {
		struct stream stream = {
			.buffer = NULL,
			.ptr = NULL,
			.in = fopen(files[i], "r")
		};
		enum efp_result res;

		if (!stream.in)
			return EFP_RESULT_FILE_NOT_FOUND;

		next_line(&stream);

		if ((res = parse_file(efp, &stream))) {
			if (stream.buffer)
				free(stream.buffer);
			fclose(stream.in);
			return res;
		}

		if (stream.buffer)
			free(stream.buffer);
		fclose(stream.in);
	}

	/* check for duplicate fragments */
	for (int i = 0; i < efp->n_lib; i++)
		for (int j = i + 1; j < efp->n_lib; j++)
			if (streq(efp->lib[i].name, efp->lib[j].name))
				return EFP_RESULT_DUPLICATE_PARAMETERS;

	for (int i = 0; i < efp->n_lib; i++) {
		set_frag_com(efp->lib + i);

		/* because of realloc we can't do this earlier */
		efp->lib[i].lib = efp->lib + i;
	}

	return EFP_RESULT_SUCCESS;
}
