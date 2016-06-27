#include "msg.h"

#ifdef WITH_MPI
#include <mpi.h>
#endif

void msg(const char *fmt, ...)
{
	va_list ap;

	va_start(ap, fmt);
	vfmsg(stdout, fmt, ap);
	va_end(ap);
}

void fmsg(FILE *st, const char *fmt, ...)
{
	va_list ap;

	va_start(ap, fmt);
	vfmsg(st, fmt, ap);
	va_end(ap);
}

void vmsg(const char *fmt, va_list ap)
{
	vfmsg(stdout, fmt, ap);
}

void vfmsg(FILE *st, const char *fmt, va_list ap)
{
#ifdef WITH_MPI
	int rank;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0) {
		vfprintf(st, fmt, ap);
	}
#else
	vfprintf(st, fmt, ap);
#endif
}
