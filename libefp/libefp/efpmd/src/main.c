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

void sim_sp(struct efp *, const struct config *);
void sim_grad(struct efp *, const struct config *);
void sim_hess(struct efp *, const struct config *);
void sim_opt(struct efp *, const struct config *);
void sim_md(struct efp *, const struct config *);

static void run_sim(struct efp *efp, const struct config *config)
{
	switch (config->run_type) {
		case RUN_TYPE_SP:   sim_sp(efp, config);   return;
		case RUN_TYPE_GRAD: sim_grad(efp, config); return;
		case RUN_TYPE_HESS: sim_hess(efp, config); return;
		case RUN_TYPE_OPT:  sim_opt(efp, config);  return;
		case RUN_TYPE_MD:   sim_md(efp, config);   return;
	}
	assert(0);
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

	return efp_strcasecmp(s1, s2);
}

static char *make_name_list(const struct config *config)
{
	size_t len = 0;

	for (int i = 0; i < config->n_frags; i++)
		len += strlen(config->frags[i].name) + 1;

	char *names = xcalloc(1, len);

	for (int i = 0; i < config->n_frags - 1; i++)
		strcat(strcat(names, config->frags[i].name), "\n");

	return strcat(names, config->frags[config->n_frags - 1].name);
}

static bool is_lib(const char *name)
{
	size_t len = strlen(name);

	return name[len - 2] == '_' && (name[len - 1] == 'l' || name[len - 1] == 'L');
}

static size_t name_len(const char *name)
{
	return is_lib(name) ? strlen(name) - 2 : strlen(name);
}

/*
 * This function constructs the list of library potential data files from
 * the list of fragment names.
 *
 * For each unique fragment: if fragment name contains an _l suffix append
 * fraglib_path prefix and remove _l suffix. Otherwise append userlib_path
 * prefix. Add ".efp" extension in both cases.
 */
static char *make_potential_file_list(const struct config *config)
{
	int i, n_uniq;
	size_t len;
	const char *uniq[config->n_frags];

	for (i = 0; i < config->n_frags; i++)
		uniq[i] = config->frags[i].name;

	qsort(uniq, config->n_frags, sizeof(char *), string_compare);

	for (i = 1, n_uniq = 1; i < config->n_frags; i++)
		if (efp_strcasecmp(uniq[i - 1], uniq[i]) != 0)
			uniq[n_uniq++] = uniq[i];

	for (i = 0, len = 0; i < n_uniq; i++) {
		const char *path = is_lib(uniq[i]) ? config->fraglib_path : config->userlib_path;
		len += strlen(path) + name_len(uniq[i]) + strlen(".efp\n") + 1;
	}

	char *files = xcalloc(1, len + 1);

	for (i = 0; i < n_uniq; i++) {
		const char *path = is_lib(uniq[i]) ? config->fraglib_path : config->userlib_path;
		strcat(strncat(strcat(strcat(files, path), "/"), uniq[i], name_len(uniq[i])), ".efp\n");
	}

	files[len - 1] = '\0';
	return files;
}

static int get_coord_count(enum efp_coord_type coord_type)
{
	switch (coord_type) {
		case EFP_COORD_TYPE_XYZABC: return  6;
		case EFP_COORD_TYPE_POINTS: return  9;
		case EFP_COORD_TYPE_ROTMAT: return 12;
	}
	assert(0);
}

static struct efp *init_sim(const struct config *config)
{
	struct efp_opts opts = {
		.terms = config->terms,
		.elec_damp = config->elec_damp,
		.disp_damp = config->disp_damp,
		.pol_damp = config->pol_damp,
		.enable_pbc = config->enable_pbc,
		.enable_cutoff = config->enable_cutoff,
		.swf_cutoff = config->swf_cutoff
	};

	char *files = make_potential_file_list(config);
	char *names = make_name_list(config);

	struct efp *efp;
	check_fail(efp_init(&efp, &opts, files, names));

	int n_coord = get_coord_count(config->coord_type);
	double coord[n_coord * config->n_frags];

	for (int i = 0; i < config->n_frags; i++)
		memcpy(coord + n_coord * i, config->frags[i].coord, n_coord * sizeof(double));

	check_fail(efp_set_coordinates(efp, config->coord_type, coord));

	if (config->enable_pbc) {
		double x = config->box[0];
		double y = config->box[1];
		double z = config->box[2];

		check_fail(efp_set_periodic_box(efp, x, y, z));
	}

	free(files);
	free(names);

	return efp;
}

static void print_usage(void)
{
	puts("usage: efpmd [-d | -v | -h | input]");
	puts("  -d  print the list of all keywords and their default values");
	puts("  -v  print package version");
	puts("  -h  print this help message");
}

int main(int argc, char **argv)
{
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

	struct config *config = parse_config(argv[1]);
	struct efp *efp = init_sim(config);

	run_sim(efp, config);

	efp_shutdown(efp);
	free_config(config);

	return EXIT_SUCCESS;
}
