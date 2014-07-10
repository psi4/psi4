import mpi4py.MPI # initialize Message Passing Interface
from ga4py import ga # initialize Global Arrays

import numpy as np

me = ga.nodeid()
nproc = ga.nnodes()

def print_distribution(g_a):
    for i in range(ga.nnodes()):
        lo,hi = ga.distribution(g_a, i)
        print "P=%s lo=%s hi=%s" % (i,lo,hi)

# create some irregular arrays
block = [3,2]
map = [0,2,6,0,5]
if nproc < np.prod(block):
    raise ValueError, "ERROR: fewer procs than requested blocks"
g_a = ga.create_irreg(ga.C_DBL, [8,10], block, map, "Array A")
if not g_a:
    ga.error("Could not create global array A",g_a)
g_b = ga.create(ga.C_INT, (2,3,4,5,6))

if not me:
    print_distribution(g_a)
    print_distribution(g_b)
