import mpi4py.MPI
from ga4py import ga
import numpy as np

me = ga.nodeid()
nproc = ga.nnodes()

def parallel_task():
    me = ga.pgroup_nodeid()
    nproc = ga.pgroup_nnodes()
    ### print a message from the master of the group
    g_a = ga.create(ga.C_DBL, (3,4,5))
    ga.randomize(g_a)
    ### sum the g_a and print the sum
    ###     -OR- do something else with g_a...

midproc = nproc//2
### assign to 'proclist_first' the first half of the process range
### assign to 'proclist_last' the last half of the process range
### create the 'group_id_first' process group
### create the 'group_id_last' process group
if me in proclist_first:
    ### set the default group to 'group_id_first'
    parallel_task()
### reset the default group to the world group
### synchronize
if me in proclist_last:
    ### set the default group to 'group_id_last'
    parallel_task()
### reset the default group to the world group
### synchronize
if not me:
    print "All done with groups"
