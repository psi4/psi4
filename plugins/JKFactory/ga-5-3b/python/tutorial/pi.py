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

### assign total number of processors to variable 'nprocs'
### assign processor ID to the variable 'myrank'

### create a global array 'g_pi' of type double and a single value

while True:
    if myrank == 0:
        n = get_n()
        ### broadcast the value of 'n'
    else:
        ### receive the broadcast of the value of 'n'
    if n == 0:
        break
    ### zero the global array 'g_pi'
    mypi = comp_pi(n, myrank, nprocs)
    ### accumulate local value 'mypi' into global array 'g_pi'
    ga.sync()
    if myrank == 0:
        ### get value of 'pi' from global array 'g_pi'
        prn_pi(pi, PI)

### destroy the global array 'g_pi'
