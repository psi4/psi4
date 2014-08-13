import mpi4py.MPI # initialize Message Passing Interface
from ga4py import ga # initialize Global Arrays

me = ga.nodeid()

def print_distribution(g_a):
    for i in range(ga.nnodes()):
        lo,hi = ga.distribution(g_a, i)
        print "%s lo=%s hi=%s" % (i,lo,hi)

# create some arrays
g_a = ga.create(ga.C_DBL, (10,20,30), chunk=(-1,20,-1))
g_b = ga.create(ga.C_DBL, (10,20,30), chunk=(10,-1,-1))

if not me:
    print_distribution(g_a)
    print_distribution(g_b)

ga.fill(g_a, 6)
ga.copy(g_a,g_b)

if not me:
    buffer = ga.access(g_b)
    print buffer.shape
    print buffer
