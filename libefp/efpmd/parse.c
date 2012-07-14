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

#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdlib.h>

#include "common.h"
#include "parse.h"

struct stream {
	char *buffer;
	char *ptr;
	FILE *in;
};

static char *read_line(FILE *in)
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
			/* fall through */
		case '\n':
			if (i == size)
				buffer = realloc(buffer, size + 1);
			buffer[i] = '\0';
			return buffer;
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

static void next_line(struct stream *stream)
{
	if (stream->buffer)
		free(stream->buffer);

	stream->buffer = read_line(stream->in);
	stream->ptr = stream->buffer;

	if (stream->buffer) {
		for (char *p = stream->buffer; *p; p++)
			*p = tolower(*p);
	}
}

static void skip_space(struct stream *stream)
{
	while (*stream->ptr && isspace(*stream->ptr))
		stream->ptr++;
}

static double read_double(struct stream *stream)
{
	skip_space(stream);

	char *endptr;
	double val = strtod(stream->ptr, &endptr);

	if (endptr == stream->ptr)
		return NAN;

	stream->ptr = endptr;
	return val;
}

static int parse_run_type(struct stream *stream, struct config *config)
{
	char *ptr = stream->ptr;

	while (*ptr && !isspace(*ptr))
		ptr++;

	config->run_type = strndup(stream->ptr, ptr - stream->ptr);
	stream->ptr = ptr;
	return 0;
}

static int parse_coord(struct stream *stream, struct config *config)
{
	static const struct {
		const char *name;
		enum coord_type value;
	} list[] = {
		{ "points", COORD_TYPE_POINTS },
		{ "xyzabc", COORD_TYPE_XYZABC }
	};

	for (size_t i = 0; i < ARRAY_SIZE(list); i++)
		if (strneq(list[i].name, stream->ptr, strlen(list[i].name))) {
			stream->ptr += strlen(list[i].name);
			config->coord_type = list[i].value;
			return 0;
		}

	return error("Unknown coord option value specified.");
}

static int parse_units(struct stream *stream, struct config *config)
{
	static const struct {
		const char *name;
		double value;
	} list[] = {
		{ "bohr", 1.0 },
		{ "angs", 1.0 / BOHR_RADIUS }
	};

	for (size_t i = 0; i < ARRAY_SIZE(list); i++)
		if (strneq(list[i].name, stream->ptr, strlen(list[i].name))) {
			stream->ptr += strlen(list[i].name);
			config->units_factor = list[i].value;
			return 0;
		}

	return error("Unknown units specified.");
}

static int parse_terms(struct stream *stream, struct config *config)
{
	static const struct {
		const char *name;
		enum efp_term value;
	} list[] = {
		{ "elec", EFP_TERM_ELEC },
		{ "pol",  EFP_TERM_POL  },
		{ "disp", EFP_TERM_DISP },
		{ "xr",   EFP_TERM_XR   }
	};

	config->efp_opts.terms = 0;

	while (*stream->ptr) {
		for (size_t i = 0; i < ARRAY_SIZE(list); i++) {
			size_t len = strlen(list[i].name);

			if (strneq(stream->ptr, list[i].name, len)) {
				stream->ptr += len;
				config->efp_opts.terms |= list[i].value;
				goto next;
			}
		}
		return error("Unknown energy term specified.");
next:
		skip_space(stream);
	}

	if (config->efp_opts.terms == 0)
		return error("At least one energy term is required.");

	return 0;
}

static int parse_disp_damp(struct stream *stream, struct config *config)
{
	static const struct {
		const char *name;
		enum efp_disp_damp value;
	} list[] = {
		{ "tt",      EFP_DISP_DAMP_TT      },
		{ "overlap", EFP_DISP_DAMP_OVERLAP },
		{ "off",     EFP_DISP_DAMP_OFF     }
	};

	for (size_t i = 0; i < ARRAY_SIZE(list); i++)
		if (strneq(list[i].name, stream->ptr, strlen(list[i].name))) {
			stream->ptr += strlen(list[i].name);
			config->efp_opts.disp_damp = list[i].value;
			return 0;
		}

	return error("Unknown dispersion damping type specified.");
}

static int parse_elec_damp(struct stream *stream, struct config *config)
{
	static const struct {
		const char *name;
		enum efp_elec_damp value;
	} list[] = {
		{ "screen",  EFP_ELEC_DAMP_SCREEN  },
		{ "overlap", EFP_ELEC_DAMP_OVERLAP },
		{ "off",     EFP_ELEC_DAMP_OFF     }
	};

	for (size_t i = 0; i < ARRAY_SIZE(list); i++)
		if (strneq(list[i].name, stream->ptr, strlen(list[i].name))) {
			stream->ptr += strlen(list[i].name);
			config->efp_opts.elec_damp = list[i].value;
			return 0;
		}

	return error("Unknown electrostatic damping type specified.");
}

static char *parse_path(struct stream *stream)
{
	char *ptr = stream->ptr;

	if (*stream->ptr == '"') {
		stream->ptr++;

		while (*stream->ptr && *stream->ptr != '"')
			stream->ptr++;

		if (!*stream->ptr)
			return NULL;
	}
	else {
		while (*stream->ptr && !isspace(*stream->ptr))
			stream->ptr++;
	}

	return strndup(ptr, stream->ptr - ptr);
}

static int parse_fraglib_path(struct stream *stream, struct config *config)
{
	free(config->fraglib_path);
	config->fraglib_path = parse_path(stream);

	if (!config->fraglib_path)
		return error("End quote not found in fraglib_path.");

	return 0;
}

static int parse_userlib_path(struct stream *stream, struct config *config)
{
	free(config->userlib_path);
	config->userlib_path = parse_path(stream);

	if (!config->userlib_path)
		return error("End quote not found in userlib_path.");

	return 0;
}

static int parse_field(struct stream *stream, struct config *config)
{
	typedef int (*parse_fn_t)(struct stream *, struct config *);

	static const struct {
		const char *name;
		parse_fn_t parse_fn;
	} list[] = {
		{ "run_type",     parse_run_type     },
		{ "coord",        parse_coord        },
		{ "units",        parse_units        },
		{ "terms",        parse_terms        },
		{ "elec_damp",    parse_elec_damp    },
		{ "disp_damp",    parse_disp_damp    },
		{ "fraglib_path", parse_fraglib_path },
		{ "userlib_path", parse_userlib_path }
	};

	for (size_t i = 0; i < ARRAY_SIZE(list); i++)
		if (strneq(list[i].name, stream->ptr, strlen(list[i].name))) {
			stream->ptr += strlen(list[i].name);
			skip_space(stream);
			return list[i].parse_fn(stream, config);
		}

	return error("Unknown parameter in input file.");
}

static int parse_frag(struct stream *stream,
		      struct config *config,
		      struct sys *sys)
{
	static const char msg[] = "Error reading fragment coordinates.";

	sys->n_frag++;
	sys->frag_name = realloc(sys->frag_name, sys->n_frag * sizeof(char *));

	skip_space(stream);

	size_t len = 0;

	while (stream->ptr[len] && !isspace(stream->ptr[len]))
		len++;

	sys->frag_name[sys->n_frag - 1] = strndup(stream->ptr, len);
	next_line(stream);

	switch (config->coord_type) {
	case COORD_TYPE_XYZABC: {
		size_t size = 6 * sys->n_frag * sizeof(double);
		sys->frag_coord = realloc(sys->frag_coord, size);
		double *ptr = sys->frag_coord + (sys->n_frag - 1) * 6;

		for (int i = 0; i < 6; i++, ptr++) {
			*ptr = read_double(stream);

			if (*ptr == NAN)
				return error(msg);

			*ptr *= config->units_factor;
		}
		next_line(stream);
		break;
	}
	case COORD_TYPE_POINTS: {
		size_t size = 9 * sys->n_frag * sizeof(double);
		sys->frag_coord = realloc(sys->frag_coord, size);
		double *ptr = sys->frag_coord + (sys->n_frag - 1) * 9;

		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++, ptr++) {
				*ptr = read_double(stream);

				if (*ptr == NAN)
					return error(msg);

				*ptr *= config->units_factor;
			}
			next_line(stream);
		}
		break;
	}
	default:
		assert(0);
	}
	return 0;
}

static int string_compare(const void *a, const void *b)
{
	const char *s1 = *(const char **)a;
	const char *s2 = *(const char **)b;

	return strcasecmp(s1, s2);
}

static char **make_potential_file_list(const char **frag_name,
				       const char *fraglib_path,
				       const char *userlib_path)
{
	/* This function constructs the list of library potential data files.
	 * For each unique fragment if fragment name contains an _l suffix
	 * append fraglib_path prefix and remove _l suffix. Otherwise append
	 * userlib_path prefix. Add .efp extension in both cases. */

	int n_frag = 0;

	while (frag_name[n_frag])
		n_frag++;

	const char **unique = malloc(n_frag * sizeof(char *));

	for (int i = 0; i < n_frag; i++)
		unique[i] = frag_name[i];

	qsort(unique, n_frag, sizeof(char *), string_compare);

	int n_unique = 1;

	for (int i = 1; i < n_frag; i++)
		if (strcasecmp(unique[i - 1], unique[i])) {
			unique[n_unique] = unique[i];
			n_unique++;
		}

	char **list = malloc((n_unique + 1) * sizeof(char *));

	for (int i = 0; i < n_unique; i++) {
		const char *name = unique[i];
		size_t len = strlen(name);
		char *path;

		if (len > 2 && name[len - 2] == '_' && name[len - 1] == 'l') {
			path = malloc(strlen(fraglib_path) + len + 4);
			strcat(strcpy(path, fraglib_path), "/");
			strcat(strncat(path, name, len - 2), ".efp");
		}
		else {
			path = malloc(strlen(userlib_path) + len + 6);
			strcat(strcpy(path, userlib_path), "/");
			strcat(strcat(path, name), ".efp");
		}

		list[i] = path;
	}

	list[n_unique] = NULL;
	free(unique);
	return list;
}

static void config_defaults(struct config *config,
			    struct sys *sys)
{
	memset(config, 0, sizeof(struct config));
	memset(sys, 0, sizeof(struct sys));

	efp_opts_default(&config->efp_opts);

	config->efp_opts.terms = EFP_TERM_ELEC | EFP_TERM_POL |
				 EFP_TERM_DISP | EFP_TERM_XR;

	config->units_factor = 1.0 / BOHR_RADIUS;
	config->fraglib_path = strdup(".");
	config->userlib_path = strdup(".");
}

int parse_config(const char *path,
		 struct config *config,
		 struct sys *sys)
{
	config_defaults(config, sys);

	struct stream stream = {
		.buffer = NULL,
		.ptr = NULL,
		.in = fopen(path, "r")
	};

	if (!stream.in)
		return error("Unable to open input file.");

	int res = 0;

	next_line(&stream);

	while (stream.buffer) {
		if (*stream.ptr == '#')
			goto next;

		skip_space(&stream);

		if (!*stream.ptr)
			goto next;

		if (strneq(stream.ptr, "fragment", 8)) {
			stream.ptr += 8;

			if ((res = parse_frag(&stream, config, sys)))
				goto fail;

			continue;
		}
		else {
			if ((res = parse_field(&stream, config)))
				goto fail;

			skip_space(&stream);

			if (*stream.ptr) {
				res = error(
			"Only one option per line is allowed.");
				goto fail;
			}
		}
next:
		next_line(&stream);
	}

	if (sys->n_frag < 1) {
		res = error("No fragments specified");
		goto fail;
	}

	sys->frag_name = realloc(sys->frag_name,
					(sys->n_frag + 1) * sizeof(char *));
	sys->frag_name[sys->n_frag] = NULL;

	config->potential_file_list =
		make_potential_file_list((const char **)sys->frag_name,
				config->fraglib_path, config->userlib_path);

	sys->gradient = malloc(6 * sys->n_frag * sizeof(double));
fail:
	fclose(stream.in);
	free(stream.buffer);
	return res;
}

void free_config(struct config *config, struct sys *sys)
{
	for (int i = 0; config->potential_file_list[i]; i++)
		free(config->potential_file_list[i]);

	free(config->potential_file_list);
	free(config->run_type);
	free(config->fraglib_path);
	free(config->userlib_path);

	for (int i = 0; sys->frag_name[i]; i++)
		free(sys->frag_name[i]);

	free(sys->frag_name);
	free(sys->frag_coord);
	free(sys->gradient);

	/* don't do free(config) and free(sys) here */
}
