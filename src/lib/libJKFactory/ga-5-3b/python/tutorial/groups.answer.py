import mpi4py.MPI
from ga4py import ga
import numpy as np

me = ga.nodeid()
nproc = ga.nnodes()

def parallel_task():
    me = ga.pgroup_nodeid()
    nproc = ga.pgroup_nnodes()
    if not me:
        print "This is process 0 on group %s" % ga.pgroup_get_default()
    g_a = ga.create(ga.C_DBL, (3,4,5))
    ga.randomize(g_a)
    if me == 0:
        print np.sum(ga.access(g_a))

midproc = nproc//2
proclist_first = range(0,midproc)
proclist_last  = range(midproc,nproc)
group_id_first = ga.pgroup_create(proclist_first)
group_id_last  = ga.pgroup_create(proclist_last)
if me in proclist_first:
    ga.pgroup_set_default(group_id_first)
    parallel_task()
ga.pgroup_set_default(ga.pgroup_get_world())
ga.sync()
if me in proclist_last:
    ga.pgroup_set_default(group_id_last)
    parallel_task()
ga.pgroup_set_default(ga.pgroup_get_world())
ga.sync()
if not me:
    print "All done with groups"
