import mpi4py.MPI # initialize Message Passing Interface
from ga4py import ga # initialize Global Arrays

print "hello from %s out of %s" % (ga.nodeid(),ga.nnodes())

