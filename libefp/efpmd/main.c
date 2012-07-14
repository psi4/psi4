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

#include <stdlib.h>

#include "common.h"
#include "parse.h"
#include "sim.h"

typedef int (*sim_fn_t)(struct efp *,
			const struct config *,
			struct sys *);

static sim_fn_t get_sim_fn(const char *run_type)
{
	static const struct {
		const char *run_type;
		sim_fn_t sim_fn;
	} sim_list[] = {
		{ "sp",   sim_sp   },
		{ "grad", sim_grad },
		{ "cg",   sim_cg   },
		{ "nve",  sim_nve  },
		{ "nvt",  sim_nvt  }
	};

	for (size_t i = 0; i < ARRAY_SIZE(sim_list); i++)
		if (streq(run_type, sim_list[i].run_type))
			return sim_list[i].sim_fn;

	return NULL;
}

static void print_banner(void)
{
	printf("EFPMD - Simple EFP simulation program based on LIBEFP\n");
	printf("Copyright (c) 2012 Ilya Kaliman\n\n");
	printf("%s\n\n", efp_banner());
}

int main(int argc, char **argv)
{
	int status = EXIT_FAILURE;

	if (argc < 2) {
		error("Usage: efpmd <input>");
		return status;
	}

	struct config config;
	struct sys sys;

	print_banner();

	if (parse_config(argv[1], &config, &sys))
		return status;

	struct efp *efp;
	enum efp_result res;
	sim_fn_t sim_fn;

	if ((res = efp_init(&efp, &config.efp_opts, NULL,
				(const char **)config.potential_file_list,
				(const char **)sys.frag_name))) {
		lib_error(res);
		goto fail;
	}

	sim_fn = get_sim_fn(config.run_type);

	if (!sim_fn) {
		error("Unknown run_type option specified.");
		goto fail;
	}

	if (sim_fn(efp, &config, &sys))
		goto fail;

	printf("EFP SIMULATION COMPLETED SUCCESSFULLY\n");
	status = EXIT_SUCCESS;
fail:
	efp_shutdown(efp);
	free_config(&config, &sys);
	return status;
}
