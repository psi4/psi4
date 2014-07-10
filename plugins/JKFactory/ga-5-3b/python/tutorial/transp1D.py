"""
transpose of 1-d array.
e.g. (1 2 3 4 5 6 7 8 9 10) => (10 9 8 7 6 5 4 3 2 1)

"""
import mpi4py.MPI
from ga4py import ga
import numpy as np

# Find local processor ID and number of processors.
### assign the local processor ID to the variable "me"
### assign the total number of processors to the variable "nprocs"

TOTALELEMS = 197
MAXPROC = 128

def verify(g_a, g_b):
    ### copy the entire block of data from the global array "g_a" into the
    ### local array "a" and similarly for "g_b" and "b".
    if not np.all(a[::-1] == b):
        print "Mismatch: a[::-1] is not equal to b"
        ga.error("verify failed")
    print "Transpose OK"

def TRANSPOSE1D():
    # Configure array dimensions. Force an unequal data distribution.
    dims = [nprocs*TOTALELEMS + nprocs/2]
    chunk = [TOTALELEMS] # minimum data on each process

    # create a global array g_a and duplicate it to get g_b
    ### create GA of integers with dimension "dims" with minimum block size
    ### "chunk" and name of "Array A" and assign the handle to the variable
    ### "g_a"
    if not g_a: ga.error("create failed: A")
    if not me: print "Created Array A"

    ### create a second global array assigned to the handled "g_b" by
    ### duplicating "g_a" and assigning the name "Array B"
    if not g_b: ga.error("duplicate failed")
    if not me: print "Created Array B"

    # initialize data in g_a
    if not me:
        print "Initializing matrix A"
        ### copy contents of a numpy range array into the remote
        ### global array "g_a"
        ### HINT: use numpy's arange() e.g. np.arange(###, dtype=np.int32)

    # Synchronize all processors to guarantee that everyone has data
    # before proceeding to the next step.
    ### synchronize all processors

    # Start initial phase of inversion by inverting the data held locally on
    # each processor. Start by finding out which data each processor owns.
    ### find out which block of data my node owns for the global array "g_a"
    ### and store the contents of the arrays into "lo" and "hi"

    # Get locally held data and copy it into local buffer a
    ### use the arrays "lo" and "hi" to copy the locally held block of data
    ### from the global array "g_a" into the local array "a".

    # Invert data locally
    b = a[::-1]

    # Invert data globally by copying locally inverted blocks into
    # their inverted positions in the GA
    lo2 = [dims[0]-hi[0]]
    hi2 = [dims[0]-lo[0]]
    ### copy data from the local array "b" into the block of the global
    ### array "g_a" described by the integer arrays "lo" and "hi"

    # Synchronize all processors to make sure inversion is complete
    ### synchronize all processors

    # Check to see if inversion is correct
    if not me: verify(g_a, g_b)

    # Deallocate arrays
    ### destroy global arrays "g_a" and "g_b"


if __name__ == '__main__':
    if not me:
        print "Using %d processes\n" % nprocs

    TRANSPOSE1D();

    if not me:
        print "\nTerminating ..."
