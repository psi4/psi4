"""Use ga.access() to sum locally per SMP node."""

import mpi4py.MPI
from ga4py import ga
import numpy as np

world_id = ga.nodeid()
world_nproc = ga.nnodes()
node_id = ga.cluster_nodeid()
node_nproc = ga.cluster_nprocs(node_id)
node_me = ga.cluster_procid(node_id,ga.nodeid())

g_a = ga.create(ga.C_DBL, (3,4,5,6))
if world_id == 0:
    ga.put(g_a, np.arange(3*4*5*6))
ga.sync()

if node_me == 0:
    sum = 0
    for i in range(node_nproc):
        smp_neighbor_world_id = ga.cluster_procid(node_id,i)
        buffer = ga.access(g_a, proc=smp_neighbor_world_id)
        sum += np.sum(buffer)
    print sum
