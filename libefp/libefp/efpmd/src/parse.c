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

#include "common.h"

#include <stream.h>

static void check_double(const struct cfg *cfg, const char *name)
{
	if (cfg_get_double(cfg, name) <= 0.0)
		error("%s must be greater than zero", name);
}

static void check_int(const struct cfg *cfg, const char *name)
{
	if (cfg_get_int(cfg, name) <= 0)
		error("%s must be greater than zero", name);
}

static void check_cfg(struct cfg *cfg)
{
	check_int(cfg, "max_steps");
	check_int(cfg, "print_step");
	check_double(cfg, "opt_tol");
	check_double(cfg, "hess_step_dist");
	check_double(cfg, "hess_step_angle");
	check_double(cfg, "time_step");
	check_double(cfg, "temperature");
	check_double(cfg, "pressure");
	check_double(cfg, "thermostat_tau");
	check_double(cfg, "barostat_tau");
}

static void parse_frag(struct stream *stream, enum efp_coord_type coord_type,
				struct frag *frag)
{
	efp_stream_skip_space(stream);
	const char *ptr = efp_stream_get_ptr(stream);
	efp_stream_skip_nonspace(stream);
	size_t len = efp_stream_get_ptr(stream) - ptr;
	frag->name = xmalloc(len + 1);

	for (size_t i = 0; i < len; i++)
		frag->name[i] = (char)tolower(*ptr++);
	frag->name[len] = '\0';

	efp_stream_next_line(stream);

	int n_rows = (int []) {
		[EFP_COORD_TYPE_XYZABC] = 1,
		[EFP_COORD_TYPE_POINTS] = 3,
		[EFP_COORD_TYPE_ROTMAT] = 4 }[coord_type];

	int n_cols = (int []) {
		[EFP_COORD_TYPE_XYZABC] = 6,
		[EFP_COORD_TYPE_POINTS] = 3,
		[EFP_COORD_TYPE_ROTMAT] = 3 }[coord_type];

	for (int i = 0, idx = 0; i < n_rows; i++) {
		for (int j = 0; j < n_cols; j++, idx++)
			if (!efp_stream_parse_double(stream, frag->coord + idx))
				error("INCORRECT FRAGMENT COORDINATES FORMAT");

		efp_stream_next_line(stream);
	}

	efp_stream_skip_space(stream);

	if (efp_stream_eol(stream))
		return;

	if (efp_strncasecmp(efp_stream_get_ptr(stream), "velocity",
			strlen("velocity")) == 0) {
		efp_stream_next_line(stream);

		for (int i = 0; i < 6; i++)
			if (!efp_stream_parse_double(stream, frag->vel + i))
				error("INCORRECT FRAGMENT VELOCITIES FORMAT");

		efp_stream_next_line(stream);
	}
}

struct sys *parse_input(struct cfg *cfg, const char *path)
{
	struct stream *stream;
	struct sys *sys = xcalloc(1, sizeof(struct sys));

	if ((stream = efp_stream_open(path)) == NULL)
		error("UNABLE TO OPEN INPUT FILE");

	efp_stream_next_line(stream);

	while (!efp_stream_eof(stream)) {
		if (efp_stream_current_char(stream) == '#')
			goto next;

		efp_stream_skip_space(stream);

		if (efp_stream_eol(stream))
			goto next;

		if (efp_strncasecmp(efp_stream_get_ptr(stream), "fragment",
				strlen("fragment")) == 0) {
			efp_stream_advance(stream, strlen("fragment"));

			sys->n_frags++;
			sys->frags = xrealloc(sys->frags, sys->n_frags * sizeof(struct frag));

			struct frag *frag = sys->frags + sys->n_frags - 1;
			memset(frag, 0, sizeof(struct frag));
			parse_frag(stream, cfg_get_enum(cfg, "coord"), frag);

			continue;
		}
		else {
			cfg_parse_line(cfg, efp_stream_get_ptr(stream));
		}
next:
		efp_stream_next_line(stream);
	}

	if (sys->n_frags < 1)
		error("AT LEAST ONE FRAGMENT MUST BE SPECIFIED");

	if (cfg_get_enum(cfg, "ensemble") == ENSEMBLE_TYPE_NPT)
		cfg_set_bool(cfg, "enable_pbc", true);

	if (cfg_get_bool(cfg, "enable_pbc"))
		cfg_set_bool(cfg, "enable_cutoff", true);

	check_cfg(cfg);
	efp_stream_close(stream);

	return sys;
}
