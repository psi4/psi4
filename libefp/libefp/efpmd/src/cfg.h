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

#ifndef EFPMD_CFG_H
#define EFPMD_CFG_H

#include <stdbool.h>
#include <efp.h>

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

struct config {
	enum run_type run_type;
	enum efp_coord_type coord_type;
	unsigned terms;
	enum efp_elec_damp elec_damp;
	enum efp_disp_damp disp_damp;
	enum efp_pol_damp pol_damp;
	bool enable_cutoff;
	double swf_cutoff;
	int max_steps;
	char *fraglib_path;
	char *userlib_path;
	bool enable_pbc;
	double box[3];
	double opt_tol;
	bool hess_central;
	double hess_step_dist;
	double hess_step_angle;
	enum ensemble_type ensemble_type;
	double time_step;
	int print_step;
	bool velocitize;
	double target_temperature;
	double target_pressure;
	double thermostat_tau;
	double barostat_tau;
	int n_frags;
	struct frag *frags;
};

struct config *parse_config(const char *);
void free_config(struct config *);
void print_defaults(void);

#endif /* EFPMD_CFG_H */
