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

#ifdef WITH_MPI
#include <mpi.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include "balance.h"
#include "private.h"

#ifdef WITH_MPI
struct master {
	int total, range[2];
};

#define MPI_CHUNK_SIZE 16

static int
master_get_work(struct master *master, int range[2])
{
#ifdef _OPENMP
#pragma omp critical
#endif
	{
		master->range[0] = master->range[1];
		master->range[1] += MPI_CHUNK_SIZE;

		if (master->range[1] > master->total)
			master->range[1] = master->total;

		range[0] = master->range[0];
		range[1] = master->range[1];
	}

	return (range[1] > range[0]);
}

static void
master_on_master(struct master *master)
{
	MPI_Status status;
	int size, range[2];

	MPI_Comm_size(MPI_COMM_WORLD, &size);

	while (master_get_work(master, range)) {
		MPI_Recv(NULL, 0, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
		MPI_Send(range, 2, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
	}

	range[0] = range[1] = -1;

	for (int i = 1; i < size; i++) {
		MPI_Recv(NULL, 0, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
		MPI_Send(range, 2, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
	}
}

static void
slave_on_master(struct master *master, struct efp *efp, work_fn fn, void *data)
{
	int range[2];

	while (master_get_work(master, range))
		fn(efp, range[0], range[1], data);
}

#ifndef _OPENMP
static int
omp_get_thread_num(void)
{
	return (0);
}
#endif /* _OPENMP */

static void
do_master(struct efp *efp, work_fn fn, void *data)
{
	struct master master;

	master.total = (int)efp->n_frag;
	master.range[0] = master.range[1] = 0;

#ifdef _OPENMP
#pragma omp parallel
#endif
	{
		if (omp_get_thread_num() == 0)
			master_on_master(&master);
		else
			slave_on_master(&master, efp, fn, data);
	}
}

static void
do_slave(struct efp *efp, work_fn fn, void *data)
{
	int range[2];

	for (;;) {
		MPI_Send(NULL, 0, MPI_INT, 0, 0, MPI_COMM_WORLD);
		MPI_Recv(range, 2, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		if (range[0] == -1 ||
		    range[1] == -1)
			break;

		fn(efp, range[0], range[1], data);
	}
}
#endif /* WITH_MPI */

void
efp_allreduce(double *x, size_t n)
{
#ifdef WITH_MPI
	MPI_Allreduce(MPI_IN_PLACE, x, (int)n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
	(void)x;
	(void)n;
#endif
}

void
efp_balance_work(struct efp *efp, work_fn fn, void *data)
{
#ifdef WITH_MPI
	int rank, size;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if (size == 1)
		fn(efp, 0, efp->n_frag, data);
	else {
		MPI_Barrier(MPI_COMM_WORLD);

		if (rank == 0)
			do_master(efp, fn, data);
		else
			do_slave(efp, fn, data);

		MPI_Barrier(MPI_COMM_WORLD);
	}
#else
	fn(efp, 0, efp->n_frag, data);
#endif
}
