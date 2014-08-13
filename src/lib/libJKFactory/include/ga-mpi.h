#ifndef GA_MPI_H_
#define GA_MPI_H_

#include <mpi.h>

extern MPI_Comm GA_MPI_Comm();
extern MPI_Comm GA_MPI_Comm_pgroup(int pgroup);
extern MPI_Comm GA_MPI_Comm_pgroup_default();

#endif /* GA_MPI_H_ */
