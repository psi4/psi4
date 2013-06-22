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

void sim_sp(struct efp *, const struct cfg *, const struct sys *);
void sim_grad(struct efp *, const struct cfg *, const struct sys *);
void sim_hess(struct efp *, const struct cfg *, const struct sys *);
void sim_opt(struct efp *, const struct cfg *, const struct sys *);
void sim_md(struct efp *, const struct cfg *, const struct sys *);

static struct cfg *make_cfg(void)
{
	struct cfg *cfg = cfg_create();

	cfg_add_enum(cfg, "run_type", RUN_TYPE_SP,
		"sp\n"
		"grad\n"
		"hess\n"
		"opt\n"
		"md\n",
		(int []) { RUN_TYPE_SP,
			   RUN_TYPE_GRAD,
			   RUN_TYPE_HESS,
			   RUN_TYPE_OPT,
			   RUN_TYPE_MD });

	cfg_add_enum(cfg, "coord", EFP_COORD_TYPE_XYZABC,
		"xyzabc\n"
		"points\n"
		"rotmat\n",
		(int []) { EFP_COORD_TYPE_XYZABC,
			   EFP_COORD_TYPE_POINTS,
			   EFP_COORD_TYPE_ROTMAT });

	cfg_add_string(cfg, "terms", "elec pol disp xr");

	cfg_add_enum(cfg, "elec_damp", EFP_ELEC_DAMP_SCREEN,
		"screen\n"
		"overlap\n"
		"off\n",
		(int []) { EFP_ELEC_DAMP_SCREEN,
			   EFP_ELEC_DAMP_OVERLAP,
			   EFP_ELEC_DAMP_OFF });

	cfg_add_enum(cfg, "disp_damp", EFP_DISP_DAMP_OVERLAP,
		"tt\n"
		"overlap\n"
		"off\n",
		(int []) { EFP_DISP_DAMP_TT,
			   EFP_DISP_DAMP_OVERLAP,
			   EFP_DISP_DAMP_OFF });

	cfg_add_enum(cfg, "pol_damp", EFP_POL_DAMP_TT,
		"tt\n"
		"off\n",
		(int []) { EFP_POL_DAMP_TT,
			   EFP_POL_DAMP_OFF });

	cfg_add_bool(cfg, "enable_cutoff", false);
	cfg_add_double(cfg, "swf_cutoff", 10.0);
	cfg_add_int(cfg, "max_steps", 100);
	cfg_add_string(cfg, "fraglib_path", FRAGLIB_PATH);
	cfg_add_string(cfg, "userlib_path", ".");
	cfg_add_bool(cfg, "enable_pbc", false);
	cfg_add_string(cfg, "periodic_box", "30.0 30.0 30.0");
	cfg_add_double(cfg, "opt_tol", 1.0e-4);
	cfg_add_bool(cfg, "hess_central", false);
	cfg_add_double(cfg, "hess_step_dist", 0.001);
	cfg_add_double(cfg, "hess_step_angle", 0.01);

	cfg_add_enum(cfg, "ensemble", ENSEMBLE_TYPE_NVE,
		"nve\n"
		"nvt\n"
		"npt\n",
		(int []) { ENSEMBLE_TYPE_NVE,
			   ENSEMBLE_TYPE_NVT,
			   ENSEMBLE_TYPE_NPT });

	cfg_add_double(cfg, "time_step", 1.0);
	cfg_add_int(cfg, "print_step", 1);
	cfg_add_bool(cfg, "velocitize", false);
	cfg_add_double(cfg, "temperature", 300.0);
	cfg_add_double(cfg, "pressure", 1.0);
	cfg_add_double(cfg, "thermostat_tau", 1.0e3);
	cfg_add_double(cfg, "barostat_tau", 1.0e4);

	return cfg;
}

static void run_sim(struct efp *efp, const struct cfg *cfg, const struct sys *sys)
{
	(void (*[])(struct efp *, const struct cfg *, const struct sys *)) {
		[RUN_TYPE_SP] = sim_sp,
		[RUN_TYPE_GRAD] = sim_grad,
		[RUN_TYPE_HESS] = sim_hess,
		[RUN_TYPE_OPT] = sim_opt,
		[RUN_TYPE_MD] = sim_md }[cfg_get_enum(cfg, "run_type")](efp, cfg, sys);
}

static void print_banner(void)
{
	printf("EFPMD ver. " LIBEFP_VERSION_STRING "\n");
	printf("Copyright (c) 2012-2013 Ilya Kaliman\n\n");
	printf("%s", efp_banner());
}

static int string_compare(const void *a, const void *b)
{
	const char *s1 = *(const char *const *)a;
	const char *s2 = *(const char *const *)b;

	return strcmp(s1, s2);
}

static bool is_lib(const char *name)
{
	size_t len = strlen(name);

	return name[len - 2] == '_' && (name[len - 1] == 'l' || name[len - 1] == 'L');
}

static void add_potentials(struct efp *efp, const struct cfg *cfg, const struct sys *sys)
{
	int i, n_uniq;
	const char *uniq[sys->n_frags];
	char path[512];

	for (i = 0; i < sys->n_frags; i++)
		uniq[i] = sys->frags[i].name;

	qsort(uniq, sys->n_frags, sizeof(const char *), string_compare);

	for (i = 1, n_uniq = 1; i < sys->n_frags; i++)
		if (strcmp(uniq[i - 1], uniq[i]) != 0)
			uniq[n_uniq++] = uniq[i];

	for (i = 0; i < n_uniq; i++) {
		const char *name = uniq[i];
		const char *prefix = is_lib(name) ?
			cfg_get_string(cfg, "fraglib_path") :
			cfg_get_string(cfg, "userlib_path");
		size_t len = is_lib(name) ? strlen(name) - 2 : strlen(name);

		snprintf(path, sizeof(path), "%s/%.*s.efp", prefix, (int)len, name);
		check_fail(efp_add_potential(efp, path));
	}
}

static unsigned get_terms(const char *str)
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

	unsigned terms = 0;

	while (*str) {
		for (size_t i = 0; i < ARRAY_SIZE(list); i++) {
			if (efp_strncasecmp(list[i].name, str, strlen(list[i].name)) == 0) {
				str += strlen(list[i].name);
				terms |= list[i].value;
				goto next;
			}
		}
		error("unknown energy term specified");
next:
		while (*str && isspace(*str))
			str++;
	}

	return terms;
}

static struct efp *init_sim(const struct cfg *cfg, const struct sys *sys)
{
	struct efp_opts opts = {
		.terms = get_terms(cfg_get_string(cfg, "terms")),
		.elec_damp = cfg_get_enum(cfg, "elec_damp"),
		.disp_damp = cfg_get_enum(cfg, "disp_damp"),
		.pol_damp = cfg_get_enum(cfg, "pol_damp"),
		.enable_pbc = cfg_get_bool(cfg, "enable_pbc"),
		.enable_cutoff = cfg_get_bool(cfg, "enable_cutoff"),
		.swf_cutoff = cfg_get_double(cfg, "swf_cutoff")
	};

	enum efp_coord_type coord = cfg_get_enum(cfg, "coord");
	struct efp *efp = efp_create();

	if (!efp)
		error("unable to create efp object");

	check_fail(efp_set_opts(efp, &opts));
	add_potentials(efp, cfg, sys);

	for (int i = 0; i < sys->n_frags; i++) {
		check_fail(efp_add_fragment(efp, sys->frags[i].name));
		check_fail(efp_set_frag_coordinates(efp, i, coord, sys->frags[i].coord));
	}

	if (cfg_get_bool(cfg, "enable_pbc")) {
		vec_t box = box_from_str(cfg_get_string(cfg, "periodic_box"));
		check_fail(efp_set_periodic_box(efp, box.x, box.y, box.z));
	}

	return efp;
}

static void print_defaults(void)
{
	struct cfg *cfg = make_cfg();
	cfg_print(cfg, stdout);
	cfg_free(cfg);
}

static void print_usage(void)
{
	puts("usage: efpmd [-d | -v | -h | input]");
	puts("  -d  print the list of all keywords and their default values");
	puts("  -v  print package version");
	puts("  -h  print this help message");
}

static void convert_units(struct cfg *cfg, struct sys *sys)
{
	cfg_set_double(cfg, "time_step",
		cfg_get_double(cfg, "time_step") * FS_TO_AU);
	cfg_set_double(cfg, "thermostat_tau",
		cfg_get_double(cfg, "thermostat_tau") * FS_TO_AU);
	cfg_set_double(cfg, "barostat_tau",
		cfg_get_double(cfg, "barostat_tau") * FS_TO_AU);
	cfg_set_double(cfg, "pressure",
		cfg_get_double(cfg, "pressure") * BAR_TO_AU);
	cfg_set_double(cfg, "swf_cutoff",
		cfg_get_double(cfg, "swf_cutoff") / BOHR_RADIUS);
	cfg_set_double(cfg, "hess_step_dist",
		cfg_get_double(cfg, "hess_step_dist") / BOHR_RADIUS);

	int n_convert = (int []) {
		[EFP_COORD_TYPE_XYZABC] = 3,
		[EFP_COORD_TYPE_POINTS] = 9,
		[EFP_COORD_TYPE_ROTMAT] = 3 }[cfg_get_enum(cfg, "coord")];

	for (int i = 0; i < sys->n_frags; i++)
		for (int j = 0; j < n_convert; j++)
			sys->frags[i].coord[j] /= BOHR_RADIUS;
}

static void sys_free(struct sys *sys)
{
	for (int i = 0; i < sys->n_frags; i++)
		free(sys->frags[i].name);

	free(sys->frags);
	free(sys);
}

int main(int argc, char **argv)
{
	struct cfg *cfg;
	struct efp *efp;
	struct sys *sys;

	if (argc < 2) {
		print_usage();
		return EXIT_FAILURE;
	}

	if (argv[1][0] == '-') {
		switch (argv[1][1]) {
		case 'd':
			print_defaults();
			return EXIT_FAILURE;
		case 'v':
			print_banner();
			return EXIT_FAILURE;
		default:
			print_usage();
			return EXIT_FAILURE;
		}
	}

	print_banner();
	printf("\n");
	cfg = make_cfg();
	sys = parse_input(cfg, argv[1]);
	printf("SIMULATION SETTINGS\n\n");
	cfg_print(cfg, stdout);
	printf("\n\n");
	convert_units(cfg, sys);
	efp = init_sim(cfg, sys);
	run_sim(efp, cfg, sys);
	efp_shutdown(efp);
	sys_free(sys);
	cfg_free(cfg);

	return EXIT_SUCCESS;
}
