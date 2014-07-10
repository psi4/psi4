"""Use ga.access() to sum locally per SMP node."""

import mpi4py.MPI
from ga4py import ga
import numpy as np

# Okay, we create the global array
g_a = ga.create(ga.C_DBL, (3,4,5,6))
if world_id == 0:
    ga.put(g_a, np.arange(3*4*5*6))
ga.sync()

# You're on your own!
