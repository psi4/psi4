from mpi4py import MPI
from ga4py import ga
from ga4py.gain import print_sync

me = ga.nodeid()
nproc = ga.nnodes()

print_sync((me,nproc))
