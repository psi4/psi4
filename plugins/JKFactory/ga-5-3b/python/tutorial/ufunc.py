import mpi4py.MPI
from ga4py import ga
import numpy as np

### I suggest creating a python "ufunc" class which implements __call__
###     but it's really up to you

g_a = ga.create(ga.C_DBL, (3,4,5))
ga.randomize(g_a)
### call a ufunc e.g. your_func(g_a)
ga.print_stdout(g_a)
### call a ufunc e.g. your_func(g_a)
ga.print_stdout(g_a)
