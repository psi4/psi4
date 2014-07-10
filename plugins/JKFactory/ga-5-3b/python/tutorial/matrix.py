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

### Assign processor ID to the int variable "me" and the total number
### of processors to the int variable "nprocs"

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
    ### create GA of doubles with dimensions "dims", with minimum block size
    ### "chunk", and with name "array A", and assign the handle to the integer
    ### variable "g_a".
    if not g_a: ga.error("create failed: A")
    if not me: print "Created Array A"

    ### Duplicate array "g_a" to create arrays "g_b" and "g_c" with array
    ### names "array B" and "array C", respectively.
    if not g_b or not g_c: ga.eror("duplicate failed")
    if not me: print "Created Arrays B and C"

    # Initialize data in matrices a and b.
    if not me: print "Initializing matrix A and B"
    a = np.random.rand(*dims)*29
    b = np.random.rand(*dims)*37

    # Copy data to global arrays g_a and g_b.
    if not me:
        ### copy the contents of array "a" into the global array "g_a"
        ### similarly for "b"

    # Synchronize all processors to make sure everyone has data.
    ### Synchronize all processors

    # Determine which block of data is locally owned. Note that
    # the same block is locally owned for all GAs.
    ### find out which block of data my node owns for the global array "g_c"
    ### and store the contents in the integer arrays "lo" and "hi"

    # Get the blocks from g_a and g_b needed to compute this block in
    # g_c and copy them into the local buffers a and b.
    lo2 = (lo[0],0)
    hi2 = (hi[0],dims[0]))
    ### copy the block of data described by the arrays "lo2" and "hi2" from
    ### the global array "g_a" in to the local array "a"

    lo3 = (0,lo[1])
    hi3 = (dims[1],hi[1]))
    ### copy the block of data described by the arrays "lo3" and "hi3" from
    ### the global array "g_b" in to the local array "b"

    # Do local matrix multiplication and store the result in local
    # buffer c. Start by evaluating the transpose of b.
    btrns = b.transpose()

    # Multiply a and b to get c.
    c = np.dot(a,b)

    # Copy c back to g_c.
    ### copy data from the local array "c" into the block of the global array
    ### "g_c" described by the integer arrays "lo" and "hi".

    verify(g_a, g_b, g_c)

    # Deallocate arrays.
    ### destroy the global arrays "g_a", "g_b", "g_c"

if __name__ == '__main__':
    if not me: print "\nUsing %d processes\n" % nprocs
    matrix_multiply()
    if not me: print "\nTerminating..."
