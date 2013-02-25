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
#include <stddef.h>

#include "common.h"

struct stream {
	char *buffer;
	char *ptr;
	FILE *in;
};

static char *read_line(FILE *in)
{
	int size = 128;
	int i = 0;
	char *buffer = xmalloc(size);

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
				buffer = xrealloc(buffer, size + 1);
			buffer[i] = '\0';
			return buffer;
		default:
			buffer[i++] = ch;

			if (i == size) {
				size *= 2;
				buffer = xrealloc(buffer, size);
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
}

static void skip_space(struct stream *stream)
{
	if (stream->ptr)
		while (*stream->ptr && isspace(*stream->ptr))
			stream->ptr++;
}

static bool parse_string(char **str, void *out)
{
	size_t len = 0;
	char *ptr = *str, *start, *outstr;

	if (!ptr)
		return false;

	while (*ptr && isspace(*ptr))
		ptr++;

	if (!*ptr)
		return false;

	if (*ptr == '"') {
		start = ++ptr;

		while (*ptr && *ptr != '"')
			ptr++, len++;

		if (!*ptr)
			return false;

		ptr++;
	}
	else {
		start = ptr;

		while (*ptr && !isspace(*ptr))
			ptr++, len++;
	}

	outstr = *(char **)out;

	if (outstr)
		free(outstr);

	outstr = xmalloc(len + 1);
	memcpy(outstr, start, len);
	outstr[len] = '\0';

	*str = ptr;
	*(char **)out = outstr;

	return true;
}

static bool parse_int(char **str, void *out)
{
	if (!*str)
		return false;

	char *endptr;
	int val = strtol(*str, &endptr, 10);

	if (endptr == *str)
		return false;

	*str = endptr;
	*(int *)out = val;
	return true;
}

static bool parse_double(char **str, void *out)
{
	if (!*str)
		return false;

	char *endptr;
	double val = strtod(*str, &endptr);

	if (endptr == *str)
		return false;

	*str = endptr;
	*(double *)out = val;
	return true;
}

static bool parse_vec(char **str, void *out)
{
	double *ptr = (double *)out;

	return parse_double(str, ptr) &&
	       parse_double(str, ptr + 1) &&
	       parse_double(str, ptr + 2);
}

static bool parse_enum(char **str, void *out, const char *names, const void *values,
				size_t data_size)
{
	const char *values_ptr = (const char *)values;

	for (const char *ptr = names; ptr; ptr = strchr(ptr, '\n')) {
		if (*ptr == '\n')
			ptr++;

		size_t len = 0;

		while (ptr[len] && ptr[len] != '\n')
			len++;

		if (efp_strncasecmp(*str, ptr, len) == 0) {
			memcpy(out, values_ptr, data_size);
			*str += len;
			return true;
		}

		values_ptr += data_size;
	}

	return false;
}

static bool parse_bool(char **str, void *out)
{
	static const char names[] =
		"true\n"
		"false";

	static const bool values[] = {
		true,
		false
	};

	return parse_enum(str, out, names, values, sizeof(values[0]));
}

static bool parse_run_type(char **str, void *out)
{
	static const char names[] =
		"sp\n"
		"grad\n"
		"hess\n"
		"opt\n"
		"md";

	static const enum run_type values[] = {
		RUN_TYPE_SP,
		RUN_TYPE_GRAD,
		RUN_TYPE_HESS,
		RUN_TYPE_OPT,
		RUN_TYPE_MD
	};

	return parse_enum(str, out, names, values, sizeof(values[0]));
}

static bool parse_coord(char **str, void *out)
{
	static const char names[] =
		"points\n"
		"xyzabc\n"
		"rotmat";

	static const enum efp_coord_type values[] = {
		EFP_COORD_TYPE_POINTS,
		EFP_COORD_TYPE_XYZABC,
		EFP_COORD_TYPE_ROTMAT
	};

	return parse_enum(str, out, names, values, sizeof(values[0]));
}

static bool parse_elec_damp(char **str, void *out)
{
	static const char names[] =
		"screen\n"
		"overlap\n"
		"off";

	static const enum efp_elec_damp values[] = {
		EFP_ELEC_DAMP_SCREEN,
		EFP_ELEC_DAMP_OVERLAP,
		EFP_ELEC_DAMP_OFF
	};

	return parse_enum(str, out, names, values, sizeof(values[0]));
}

static bool parse_disp_damp(char **str, void *out)
{
	static const char names[] =
		"tt\n"
		"overlap\n"
		"off";

	static const enum efp_disp_damp values[] = {
		EFP_DISP_DAMP_TT,
		EFP_DISP_DAMP_OVERLAP,
		EFP_DISP_DAMP_OFF
	};

	return parse_enum(str, out, names, values, sizeof(values[0]));
}

static bool parse_pol_damp(char **str, void *out)
{
	static const char names[] =
		"tt\n"
		"off";

	static const enum efp_pol_damp values[] = {
		EFP_POL_DAMP_TT,
		EFP_POL_DAMP_OFF
	};

	return parse_enum(str, out, names, values, sizeof(values[0]));
}

static bool parse_ensemble(char **str, void *out)
{
	static const char names[] =
		"nve\n"
		"nvt\n"
		"npt";

	static const enum ensemble_type values[] = {
		ENSEMBLE_TYPE_NVE,
		ENSEMBLE_TYPE_NVT,
		ENSEMBLE_TYPE_NPT
	};

	return parse_enum(str, out, names, values, sizeof(values[0]));
}

static bool parse_terms(char **str, void *out)
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

	char *ptr = *str;
	unsigned terms = 0;

	while (*ptr) {
		for (size_t i = 0; i < ARRAY_SIZE(list); i++) {
			if (efp_strncasecmp(list[i].name, ptr, strlen(list[i].name)) == 0) {
				ptr += strlen(list[i].name);
				terms |= list[i].value;
				goto next;
			}
		}
		return false;
next:
		while (*ptr && isspace(*ptr))
			ptr++;
	}

	*(unsigned *)out = terms;
	*str = ptr;
	return terms != 0;
}

static bool int_gt_zero(void *val)
{
	return *(int *)val > 0;
}

static bool double_gt_zero(void *val)
{
	return *(double *)val > 0.0;
}

static const struct {
	const char *name;
	const char *default_value;
	bool (*parse_fn)(char **, void *);
	bool (*check_fn)(void *);
	size_t member_offset;
} config_list[] = {
#include "keywords.h"
};

static void convert_units(struct config *config)
{
	config->time_step *= FS_TO_AU;
	config->thermostat_tau *= FS_TO_AU;
	config->barostat_tau *= FS_TO_AU;
	config->target_pressure *= BAR_TO_AU;
	config->box[0] /= BOHR_RADIUS;
	config->box[1] /= BOHR_RADIUS;
	config->box[2] /= BOHR_RADIUS;
	config->swf_cutoff /= BOHR_RADIUS;
	config->hess_step_dist /= BOHR_RADIUS;

	int n_convert = 0;

	switch (config->coord_type) {
		case EFP_COORD_TYPE_XYZABC: n_convert = 3; break;
		case EFP_COORD_TYPE_POINTS: n_convert = 9; break;
		case EFP_COORD_TYPE_ROTMAT: n_convert = 3; break;
	};

	for (int i = 0; i < config->n_frags; i++)
		for (int j = 0; j < n_convert; j++)
			config->frags[i].coord[j] /= BOHR_RADIUS;
}

static void parse_field(struct stream *stream, struct config *config)
{
	for (size_t i = 0; i < ARRAY_SIZE(config_list); i++) {
		const char *name = config_list[i].name;
		size_t offset = config_list[i].member_offset;

		if (efp_strncasecmp(name, stream->ptr, strlen(name)) == 0) {
			if (!isspace(stream->ptr[strlen(name)]))
				continue;

			stream->ptr += strlen(name);
			skip_space(stream);

			if (!config_list[i].parse_fn(&stream->ptr, (char *)config + offset))
				error("INCORRECT VALUE FOR KEYWORD '%s'", name);

			bool (*check_fn)(void *) = config_list[i].check_fn;

			if (check_fn && !check_fn((char *)config + offset))
				error("VALUE IS OUT OF RANGE FOR KEYWORD '%s'", name);

			return;
		}
	}

	int len = 0;

	while (stream->ptr[len] && !isspace(stream->ptr[len]))
		len++;

	error("UNKNOWN KEYWORD '%.*s'", len, stream->ptr);
}

static void parse_frag(struct stream *stream, enum efp_coord_type coord_type,
				struct frag *frag)
{
	memset(frag, 0, sizeof(struct frag));

	if (!parse_string(&stream->ptr, &frag->name))
		error("UNABLE TO READ FRAGMENT NAME");

	for (char *p = frag->name; *p; p++)
		*p = tolower(*p);

	next_line(stream);

	int n_rows = 0, n_cols = 0;

	switch (coord_type) {
		case EFP_COORD_TYPE_XYZABC: n_rows = 1; n_cols = 6; break;
		case EFP_COORD_TYPE_POINTS: n_rows = 3; n_cols = 3; break;
		case EFP_COORD_TYPE_ROTMAT: n_rows = 4; n_cols = 3; break;
	}

	for (int i = 0, idx = 0; i < n_rows; i++) {
		for (int j = 0; j < n_cols; j++, idx++)
			if (!parse_double(&stream->ptr, frag->coord + idx))
				error("INCORRECT FRAGMENT COORDINATES FORMAT");

		next_line(stream);
	}

	if (!stream->ptr)
		return;

	skip_space(stream);

	if (efp_strncasecmp(stream->ptr, "velocity", strlen("velocity")) == 0) {
		next_line(stream);

		for (int i = 0; i < 6; i++)
			if (!parse_double(&stream->ptr, frag->vel + i))
				error("INCORRECT FRAGMENT VELOCITIES FORMAT");

		next_line(stream);
	}
}

static void set_config_defaults(struct config *config)
{
	memset(config, 0, sizeof(struct config));

	for (size_t i = 0; i < ARRAY_SIZE(config_list); i++) {
		size_t offset = config_list[i].member_offset;

		char buffer[128], *ptr = buffer;
		strncpy(buffer, config_list[i].default_value, sizeof(buffer));

		bool res = config_list[i].parse_fn(&ptr, (char *)config + offset);
		assert(res);
	}
}

struct config *parse_config(const char *path)
{
	struct config *config = xmalloc(sizeof(struct config));
	set_config_defaults(config);

	struct stream stream = {
		.buffer = NULL,
		.ptr = NULL,
		.in = fopen(path, "r")
	};

	if (!stream.in)
		error("UNABLE TO OPEN INPUT FILE");

	next_line(&stream);

	while (stream.buffer) {
		if (*stream.ptr == '#')
			goto next;

		skip_space(&stream);

		if (!*stream.ptr)
			goto next;

		if (efp_strncasecmp(stream.ptr, "fragment", strlen("fragment")) == 0) {
			stream.ptr += strlen("fragment");

			config->n_frags++;
			config->frags = xrealloc(config->frags,
				config->n_frags * sizeof(struct frag));

			struct frag *frag = config->frags + config->n_frags - 1;
			parse_frag(&stream, config->coord_type, frag);

			continue;
		}
		else {
			parse_field(&stream, config);
			skip_space(&stream);

			if (*stream.ptr)
				error("ONLY ONE OPTION PER LINE IS ALLOWED");
		}
next:
		next_line(&stream);
	}

	if (config->n_frags < 1)
		error("AT LEAST ONE FRAGMENT MUST BE SPECIFIED");

	if (config->ensemble_type == ENSEMBLE_TYPE_NPT)
		config->enable_pbc = true;

	if (config->enable_pbc)
		config->enable_cutoff = true;

	convert_units(config);

	fclose(stream.in);
	free(stream.buffer);

	return config;
}

void free_config(struct config *config)
{
	free(config->fraglib_path);
	free(config->userlib_path);

	for (int i = 0; i < config->n_frags; i++)
		free(config->frags[i].name);

	free(config->frags);
	free(config);
}

void print_defaults(void)
{
	for (size_t i = 0; i < ARRAY_SIZE(config_list); i++)
		printf("%s %s\n", config_list[i].name, config_list[i].default_value);
}
