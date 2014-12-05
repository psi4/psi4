/*-
 * Copyright (c) 2012-2014 Ilya Kaliman
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

#include <stdarg.h>
#include <stdio.h>

#ifdef WITH_MPI
#include <mpi.h>
#endif

#include "log.h"

static efp_log_cb _log_cb = NULL;

void
efp_log(const char *fmt, ...)
{
	va_list ap;
	char msg[512];

	if (_log_cb == NULL)
		return;

#ifdef _OPENMP
#pragma omp master
#endif
{
#ifdef WITH_MPI
	int rank;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0) {
		va_start(ap, fmt);
		vsnprintf(msg, sizeof(msg), fmt, ap);
		_log_cb(msg);
		va_end(ap);
	}
#else
	va_start(ap, fmt);
	vsnprintf(msg, sizeof(msg), fmt, ap);
	_log_cb(msg);
	va_end(ap);
#endif
}
}

void
efp_set_log_cb(efp_log_cb log_cb)
{
	_log_cb = log_cb;
}

efp_log_cb
efp_get_log_cb(void)
{
	return (_log_cb);
}
