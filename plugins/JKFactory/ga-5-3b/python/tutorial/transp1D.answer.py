"""
transpose of 1-d array.
e.g. (1 2 3 4 5 6 7 8 9 10) => (10 9 8 7 6 5 4 3 2 1)

"""
import mpi4py.MPI
from ga4py import ga
import numpy as np

# Find local processor ID and number of processors.
me = ga.nodeid()
nprocs = ga.nnodes()

TOTALELEMS = 197

def verify(g_a, g_b):
    a = ga.get(g_a)
    b = ga.get(g_b)
    if not np.all(a[::-1] == b):
        print "Mismatch: a[::-1] is not equal to b"
        ga.error("verify failed")
    print "Transpose OK"

def TRANSPOSE1D():
    # Configure array dimensions. Force an unequal data distribution.
    dims = [nprocs*TOTALELEMS + nprocs/2]
    chunk = [TOTALELEMS] # minimum data on each process

    # create a global array g_a and duplicate it to get g_b
    g_a = ga.create(ga.C_INT, dims, "array A", chunk)
    if not g_a: ga.error("create failed: A")
    if not me: print "Created Array A"

    g_b = ga.duplicate(g_a, "array B")
    if not g_b: ga.error("duplicate failed")
    if not me: print "Created Array B"

    # initialize data in g_a
    if not me:
        print "Initializing matrix A"
        ga.put(g_a, np.arange(dims[0], dtype=np.int32))

    # Synchronize all processors to guarantee that everyone has data
    # before proceeding to the next step.
    ga.sync()

    # Start initial phase of inversion by inverting the data held locally on
    # each processor. Start by finding out which data each processor owns.
    lo,hi = ga.distribution(g_a)

    # Get locally held data and copy it into local buffer a
    a = ga.get(g_a, lo, hi)

    # Invert data locally
    b = a[::-1]

    # Invert data globally by copying locally inverted blocks into
    # their inverted positions in the GA
    ga.put(g_b, b, dims[0]-hi[0], dims[0]-lo[0])

    # Synchronize all processors to make sure inversion is complete
    ga.sync()

    # Check to see if inversion is correct
    if not me: verify(g_a, g_b)

    # Deallocate arrays
    ga.destroy(g_a)
    ga.destroy(g_b)


if __name__ == '__main__':
    if not me:
        print "Using %d processes\n" % nprocs

    TRANSPOSE1D();

    if not me:
        print "\nTerminating ..."
