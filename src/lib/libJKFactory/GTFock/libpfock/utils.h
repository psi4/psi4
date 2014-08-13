#ifndef __UTILS_H__
#define __UTILS_H__


void printmatRM (char *name, double *A, int nrows, int ncols);

void printmatCM (char *name, double *A, int nrows, int ncols);

void initomp (int nthreads, int verbose);

void *my_malloc (size_t size, char *file, int line);


#endif /* __UTILS_H__ */
