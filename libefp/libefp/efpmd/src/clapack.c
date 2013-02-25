#include <stdlib.h>

#include "clapack.h"

void dsyev_(char *,
	    char *,
	    int *,
	    double *,
	    int *,
	    double *,
	    double *,
	    int *,
	    int *);

int c_dsyev(char jobz, char uplo, int n, double *a, int lda, double *w)
{
	int info, lwork;
	double *work;

	lwork = n * n;
	work = malloc(lwork * sizeof(double));

	dsyev_(&jobz, &uplo, &n, a, &lda, w, work, &lwork, &info);

	free(work);
	return info;
}
