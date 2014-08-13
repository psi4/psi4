import mpi4py.MPI
from ga4py import ga
import numpy as np

class ufunc(object):
    def __init__(self, ufunc):
        self.ufunc = ufunc
    def __call__(self, g_a):
        a = ga.access(g_a)
        if a is not None:
            self.ufunc(a,a)

sin = ufunc(np.sin)
cos = ufunc(np.cos)

g_a = ga.create(ga.C_DBL, (3,4,5))
ga.randomize(g_a)
sin(g_a)
ga.print_stdout(g_a)
cos(g_a)
ga.print_stdout(g_a)
