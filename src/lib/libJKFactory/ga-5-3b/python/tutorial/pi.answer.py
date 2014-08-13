#!/usr/bin/env python
"""
THIS CODE WAS TAKEN FROM mpi4py-1.2.2/demo/compute-pi and converted to use
Global Arrays and ga.acc().  Thanks Lisandro!

Parallel PI computation using Remote Memory Access (RMA)
within Python objects exposing memory buffers (requires NumPy).

usage::

  $ mpiexec -n <nprocs> python pi.py

"""
from mpi4py import MPI
from ga4py import ga

from math   import pi as PI
from numpy  import ndarray
import sys

def get_n():
    prompt  = "Enter the number of intervals: (0 quits) "
    try:
        n = int(raw_input(prompt));
        if n < 0: n = 0
    except:
        n = 0
    return n

def comp_pi(n, myrank=0, nprocs=1):
    h = 1.0 / n;
    s = 0.0;
    for i in xrange(myrank + 1, n + 1, nprocs):
        x = h * (i - 0.5);
        s += 4.0 / (1.0 + x**2);
    return s * h

def prn_pi(pi, PI):
    message = "pi is approximately %.16f, error is %.16f"
    print  (message % (pi, abs(pi - PI)))

nprocs = ga.nnodes()
myrank = ga.nodeid()

g_pi = ga.create(ga.C_DBL, [1])

one_time = False
if len(sys.argv) == 2:
    n = int(sys.argv[1])
    one_time = True

while True:
    if not one_time:
        if myrank == 0:
            n = get_n()
            n = ga.brdcst(n)
        else:
            n = ga.brdcst(0)
        if n == 0:
            break
    ga.zero(g_pi)
    mypi = comp_pi(n, myrank, nprocs)
    ga.acc(g_pi, mypi)
    ga.sync()
    if myrank == 0:
        pi = ga.get(g_pi)[0]
        prn_pi(pi, PI)
    if one_time:
        break

ga.destroy(g_pi)
