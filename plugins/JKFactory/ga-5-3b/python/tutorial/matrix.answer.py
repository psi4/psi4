"""
Multiplication of two square matrices with randomly generated contents.

"""
import mpi4py.MPI
from ga4py import ga

import numpy as np

NDIM = 2
TOTALELEMS = 1007
MAXPROC = 128
NBUF = 4
TOLERANCE = 0.1

me = ga.nodeid()
nprocs = ga.nnodes()

def verify(g_a, g_b, g_c):
    g_chk = ga.duplicate(g_a, "array check")
    if not g_chk: ga.error("duplicate failed")
    ga.sync()

    ga.gemm(False, False, TOTALELEMS, TOTALELEMS, TOTALELEMS, 1.0, g_a, g_b,
            0.0, g_chk);
    ga.sync()

    ga.add(g_c, g_chk, g_chk, 1.0, -1.0)
    rchk = ga.dot(g_chk, g_chk)

    if not me:
        print "Normed difference in matrices: %12.4f" % rchk
        if not (-TOLERANCE < rchk < TOLERANCE):
            ga.error("Matrix multiply verify failed")
        else:
            print "Matrix Multiply OK"

    ga.destroy(g_chk)

def matrix_multiply():
    # Configure array dimensions. Force an unequal data distribution.
    dims = [TOTALELEMS]*NDIM
    chunk = [TOTALELEMS/nprocs-1]*NDIM

    # Create a global array g_a and duplicate it to get g_b and g_c.
    g_a = ga.create(ga.C_DBL, dims, "array A", chunk)
    if not g_a: ga.error("create failed: A")
    if not me: print "Created Array A"

    g_b = ga.duplicate(g_a, "array B")
    g_c = ga.duplicate(g_a, "array C")
    if not g_b or not g_c: ga.eror("duplicate failed")
    if not me: print "Created Arrays B and C"

    # Initialize data in matrices a and b.
    if not me: print "Initializing matrix A and B"
    a = np.random.rand(*dims)*29
    b = np.random.rand(*dims)*37

    # Copy data to global arrays g_a and g_b.
    if not me:
        ga.put(g_a, a)
        ga.put(g_b, b)

    # Synchronize all processors to make sure everyone has data.
    ga.sync()

    # Determine which block of data is locally owned. Note that
    # the same block is locally owned for all GAs.
    lo,hi = ga.distribution(g_c)

    # Get the blocks from g_a and g_b needed to compute this block in
    # g_c and copy them into the local buffers a and b.
    a = ga.get(g_a, (lo[0],0), (hi[0],dims[0]))
    b = ga.get(g_b, (0,lo[1]), (dims[1],hi[1]))

    # Do local matrix multiplication and store the result in local
    # buffer c. Start by evaluating the transpose of b.
    btrns = b.transpose()

    # Multiply a and b to get c.
    c = np.dot(a,b)

    # Copy c back to g_c.
    ga.put(g_c, c, lo, hi)

    verify(g_a, g_b, g_c)

    # Deallocate arrays.
    ga.destroy(g_a)
    ga.destroy(g_b)
    ga.destroy(g_c)

if __name__ == '__main__':
    if not me: print "\nUsing %d processes\n" % nprocs
    matrix_multiply()
    if not me: print "\nTerminating..."
