#!/usr/bin/env python
"""At first it was a duplicate of GA's test.x, but was made more modular since
so much could be reused.  It is more like unit tests since each piece has
build up and tear down."""

import random
import sys

from mpi4py import MPI
from ga4py import ga
import numpy as np

me = ga.nodeid()
nproc = ga.nnodes()
nnodes = ga.cluster_nnodes()
inode = ga.cluster_nodeid()
lprocs = ga.cluster_nprocs(inode)
iproc = me % lprocs

n = 256
m = 2*n
maxloop = 100
nloop = min(maxloop,n)
block_size = [32,32]
proc_grid = [2,nproc/2]
MEM_INC = 1000
MIRROR = False
USE_RESTRICTED = False
NEW_API = False
NGA_GATSCAT = False
BLOCK_CYCLIC = False
USE_SCALAPACK_DISTR = False

# This was used to debug as it makes printed numpy arrays a bit easier to read
np.set_printoptions(precision=6, suppress=True, edgeitems=4)

def mismatch(x,y):
    #return abs(x-y)/max(1.0,abs(x)) > 1e-12
    return abs(x-y)/max(1.0,abs(x)) > 1e-5

def main():
    if 0 == me:
        if MIRROR:
            print ' Performing tests on Mirrored Arrays'
        print ' GA initialized'

    # note that MA is not used, so no need to initialize it
    # "import ga" registers malloc/free as memory allocators internally

    #if nproc-1 == me:
    if 0 == me:
        print 'using %d process(es) %d custer nodes' % (
                nproc, ga.cluster_nnodes())
        print 'process %d is on node %d with %d processes' % (
                me, ga.cluster_nodeid(), ga.cluster_nprocs(-1))

    # create array to force staggering of memory and uneven distribution
    # of pointers
    dim1 = MEM_INC
    mapc = [0]*nproc
    for i in range(nproc):
        mapc[i] = MEM_INC*i
        dim1 += MEM_INC*i
    g_s = ga.create_handle()
    ga.set_data(g_s, [dim1], ga.C_INT)
    ga.set_array_name(g_s, 's')
    ga.set_irreg_distr(g_s, mapc, [nproc])

    if MIRROR:
        if 0 == me:
            print '\nTESTING MIRRORED ARRAYS\n'

    # check support for single precision arrays
    if 0 == me:
        print '\nCHECKING SINGLE PRECISION\n'
    check_float()

    # check support for double precision arrays
    if 0 == me:
        print '\nCHECKING DOUBLE PRECISION\n'
    check_double()

    # check support for single precision complex arrays
    if 0 == me:
        print '\nCHECKING SINGLE COMPLEX\n'
    check_complex_float()

    # check support for double precision complex arrays
    if 0 == me:
        print '\nCHECKING DOUBLE COMPLEX\n'
    check_complex_double()

    # check support for integer arrays
    if 0 == me:
        print '\nCHECKING INT\n'
    check_int()

    # check support for long integer arrays
    if 0 == me:
        print '\nCHECKING LONG INT\n'
    check_long()

    if 0 == me:
        print '\nCHECKING Wrappers to Message Passing Collective ops\n'
    check_wrappers()

    # check if memory limits are enforced
    #check_mem(ma_heap*ga.nnodes())

    if 0 == me: ga.print_stats()
    if 0 == me: print ' '
    if 0 == me: print 'All tests successful'

    # tidy up the ga package
    # NO NEED -- atexit registered ga.terminate()

    # tidy up after message-passing library
    # NO NEED -- atexit registered MPI.Finalize()

    # Note: so long as mpi4py is imported before ga, cleanup is automatic

def create_local_a(gatype):
    nptype = ga.dtype(gatype)
    if gatype == ga.C_SCPL:
        if MIRROR:
            a = np.fromfunction(lambda i,j: i+inode, (n,n), dtype=np.float32)
            b = np.fromfunction(lambda i,j: j*n,     (n,n), dtype=np.float32)
            return np.vectorize(complex)(a,b)
        else:
            a = np.fromfunction(lambda i,j: i,   (n,n), dtype=np.float32)
            b = np.fromfunction(lambda i,j: j*n, (n,n), dtype=np.float32)
            return np.vectorize(complex)(a,b)
    elif gatype == ga.C_DCPL:
        if MIRROR:
            a = np.fromfunction(lambda i,j: i+inode, (n,n), dtype=np.float64)
            b = np.fromfunction(lambda i,j: j*n,     (n,n), dtype=np.float64)
            return np.vectorize(complex)(a,b)
        else:
            a = np.fromfunction(lambda i,j: i,   (n,n), dtype=np.float64)
            b = np.fromfunction(lambda i,j: j*n, (n,n), dtype=np.float64)
            return np.vectorize(complex)(a,b)
    elif gatype in [ga.C_DBL,ga.C_FLOAT]:
        if MIRROR:
            return np.fromfunction(lambda i,j: inode+i+j*n, (n,n), dtype=nptype)
        else:
            return np.fromfunction(lambda i,j: i+j*n,       (n,n), dtype=nptype)
    elif gatype in [ga.C_INT,ga.C_LONG]:
        if MIRROR:
            return np.fromfunction(lambda i,j: inode+i+j*1000, (n,n), dtype=nptype)
        else:
            return np.fromfunction(lambda i,j: i+j*1000, (n,n), dtype=nptype)

def create_local_b(gatype):
    nptype = ga.dtype(gatype)
    b = np.zeros((n,n), dtype=nptype)
    if gatype in [ga.C_SCPL, ga.C_DCPL]:
        b[:] = complex(-1,1)
    else:
        b[:] = -1
    return b

def create_global_array(gatype):
    if NEW_API:
        g_a = ga.create_handle()
        ga.set_data(g_a, [n,n], gatype)
        ga.set_array_name(g_a, 'a')
        if USE_RESTRICTED:
            num_restricted = nproc/2 or 1
            restricted_list = np.arange(num_restricted) + num_restricted/2
            ga.set_restricted(g_a, restricted_list)
        if BLOCK_CYCLIC:
            if USE_SCALAPACK_DISTR:
                if nproc % 2 == 0:
                    ga.error('Available procs must be divisible by 2',nproc)
                ga.set_block_cyclic_proc_grid(g_a, block_size, proc_grid)
            else:
                ga.set_block_cyclic(g_a, block_size)
        if MIRROR:
            p_mirror = ga.pgroup_get_mirror()
            ga.set_pgroup(g_a, p_mirror)
        ga.allocate(g_a)
    else:
        if MIRROR:
            p_mirror = ga.pgroup_get_mirror()
            ga.create_config(gatype, (n,n), 'a', None, p_mirror)
        else:
            g_a = ga.create(gatype, (n,n), 'a')
    if 0 == g_a:
        ga.error('ga.create failed')
    if MIRROR:
        lproc = me - ga.cluster_procid(inode, 0)
        lo,hi = ga.distribution(g_a, lproc)
    else:
        lo,hi = ga.distribution(g_a, me)
    ga.sync()
    return g_a

def check_zero(gatype):
    if 0 == me:
        print '> Checking zero ...',
    g_a = create_global_array(gatype)
    ga.zero(g_a)
    a = ga.get(g_a)
    if not np.all(a == 0):
        ga.error('ga.zero failed')
    if 0 == me:
        print 'OK'
    ga.destroy(g_a)

def check_put_disjoint(gatype):
    """each node fills in disjoint sections of the array"""
    if 0 == me:
        print '> Checking disjoint put ...',
    g_a = create_global_array(gatype)
    a = create_local_a(gatype)
    inc = (n-1)/20 + 1
    ij = 0
    for i in range(0,n,inc):
        for j in range(0,n,inc):
            check = False
            if MIRROR:
                check = ij % lprocs == iproc
            else:
                check = ij % nproc == me
            if check:
                lo = [i,j]
                hi = [min(i+inc,n), min(j+inc,n)]
                piece = a[ga.zip(lo,hi)]
                ga.put(g_a, piece, lo, hi)
                # the following check is not part of the original test.F
                result = ga.get(g_a, lo, hi)
                if not np.all(result == piece):
                    ga.error("put followed by get failed", 1)
            ga.sync()
            ij += 1
    ga.sync()
    # all nodes check all of a
    b = ga.get(g_a)
    if not np.all(a == b):
        ga.error('put failed, exiting')
    if 0 == me:
        print 'OK'
    ga.destroy(g_a)

def check_get(gatype):
    """check nloop random gets from each node"""
    if 0 == me:
        print '> Checking random get (%d calls)...' % nloop
    g_a = create_global_array(gatype)
    a = create_local_a(gatype)
    if 0 == me:
        ga.put(g_a, a)
    ga.sync()
    nwords = 0
    random.seed(ga.nodeid()*51+1) # different seed for each proc
    for loop in range(nloop):
        ilo,ihi = random.randint(0, nloop-1),random.randint(0, nloop-1)
        if ihi < ilo: ilo,ihi = ihi,ilo
        jlo,jhi = random.randint(0, nloop-1),random.randint(0, nloop-1)
        if jhi < jlo: jlo,jhi = jhi,jlo
        nwords += (ihi-ilo+1)*(jhi-jlo+1)
        ihi += 1
        jhi += 1
        result = ga.get(g_a, (ilo,jlo), (ihi,jhi))
        if not np.all(result == a[ilo:ihi,jlo:jhi]):
            ga.error('random get failed')
        if 0 == me and loop % max(1,nloop/20) == 0:
            print ' call %d node %d checking get((%d,%d),(%d,%d)) total %f' % (
                    loop, me, ilo, ihi, jlo, jhi, nwords)
    if 0 == me:
        print 'OK'
    ga.destroy(g_a)

def check_accumulate_disjoint(gatype):
    """Each node accumulates into disjoint sections of the array."""
    if 0 == me:
        print '> Checking disjoint accumulate ...',
    g_a = create_global_array(gatype)
    a = create_local_a(gatype)
    b = np.fromfunction(lambda i,j: i+j+2, (n,n), dtype=ga.dtype(gatype))
    if 0 == me:
        ga.put(g_a, a)
    ga.sync()
    inc = (n-1)/20 + 1
    ij = 0
    for i in range(0,n,inc):
        for j in range(0,n,inc):
            x = 10.0
            lo = [i,j]
            hi = [min(i+inc,n), min(j+inc,n)]
            piece = b[ga.zip(lo,hi)]
            check = False
            if MIRROR:
                check = ij % lprocs == iproc
            else:
                check = ij % nproc == me
            if check:
                ga.acc(g_a, piece, lo, hi, x)
            ga.sync()
            ij += 1
            # each process applies all updates to its local copy
            a[ga.zip(lo,hi)] += x * piece
    ga.sync()
    # all nodes check all of a
    if not np.all(ga.get(g_a) == a):
        ga.error('acc failed')
    if 0 == me:
        print 'OK'
    ga.destroy(g_a)

def check_accumulate_overlap(gatype):
    if 0 == me:
        print '> Checking overlapping accumulate ...',
    g_a = create_global_array(gatype)
    ga.zero(g_a)
    ga.acc(g_a, [1], (n/2,n/2), (n/2+1,n/2+1), 1)
    ga.sync()
    if MIRROR:
        if 0 == iproc:
            x = abs(ga.get(g_a, (n/2,n/2), (n/2+1,n/2+1))[0,0] - lprocs)
            if not 0 == x:
                ga.error('overlapping accumulate failed -- expected %s got %s'%(
                        x, lprocs))
    else:
        if 0 == me:
            x = abs(ga.get(g_a, (n/2,n/2), (n/2+1,n/2+1))[0,0] - nproc)
            if not 0 == x:
                ga.error('overlapping accumulate failed -- expected %s got %s'%(
                        x, nproc))
    if 0 == me:
        print 'OK'
    ga.destroy(g_a)

def check_add(gatype):
    if 0 == me:
        print '> Checking add ...',
    g_a = create_global_array(gatype)
    g_b = create_global_array(gatype)
    a = create_local_a(gatype)
    b = create_local_b(gatype)
    alpha = None
    beta = None
    if 0 == me:
        ga.put(g_a, a)
    ga.sync();
    np.random.seed(12345) # everyone has same seed
    if gatype in [ga.C_SCPL,ga.C_DCPL]:
        b_real = np.random.random_sample((n,n))
        b_imag = np.random.random_sample((n,n))
        b[:] = np.vectorize(complex)(b_real,b_imag)
        alpha = complex(0.1,-0.1)
        beta = complex(0.9,-0.9)
    else:
        b[:] = np.random.random_sample((n,n))
        alpha = 0.1
        beta = 0.9
    a = alpha*a + beta*b
    if MIRROR:
        if 0 == iproc:
            ga.put(g_b, b)
    else:
        if 0 == me:
            ga.put(g_b, b)
    ga.sync()
    ga.add(g_a, g_b, g_b, alpha, beta)
    b = ga.get(g_b, buffer=b)
    if np.any(np.vectorize(mismatch)(b,a)):
        ga.error('add failed')
    if 0 == me:
        print 'OK'
    ga.destroy(g_a)
    ga.destroy(g_b)

def check_dot(gatype):
    if 0 == me:
        print '> Checking dot ...',
    np.random.seed(12345) # everyone has same seed
    g_a = create_global_array(gatype)
    g_b = create_global_array(gatype)
    a = create_local_a(gatype)
    b = np.random.random_sample((n,n))
    if MIRROR:
        if 0 == iproc:
            ga.put(g_b, b)
            ga.put(g_a, a)
    else:
        if 0 == me:
            ga.put(g_b, b)
            ga.put(g_a, a)
    ga.sync()
    sum1 = np.sum(a*b)
    sum2 = ga.dot(g_a, g_b)
    if mismatch(sum1, sum2):
        ga.error('dot wrong %s != %s' % (sum1, sum2))
    if 0 == me:
        print 'OK'
    ga.destroy(g_a)
    ga.destroy(g_b)

def check_scale(gatype):
    if 0 == me:
        print '> Checking scale ...',
    g_a = create_global_array(gatype)
    a = create_local_a(gatype)
    if 0 == me:
        ga.put(g_a, a)
    ga.sync()
    ga.scale(g_a, 0.123)
    a *= 0.123
    if np.any(np.vectorize(mismatch)(a,ga.get(g_a))):
        ga.error('add failed')
    if 0 == me:
        print 'OK'
    ga.destroy(g_a)

def check_copy(gatype):
    if 0 == me:
        print '> Checking copy ...',
    g_a = create_global_array(gatype)
    g_b = create_global_array(gatype)
    a = create_local_a(gatype)
    if 0 == me:
        ga.put(g_a, a)
    ga.copy(g_a, g_b)
    if not np.all(a == ga.get(g_b)):
        ga.error('copy failed')
    if 0 == me:
        print 'OK'
    ga.destroy(g_a)
    ga.destroy(g_b)

def check_gather(gatype):
    if 0 == me:
        print '> Checking gather (might be slow)...',
    g_a = create_global_array(gatype)
    a = create_local_a(gatype)
    if 0 == me:
        ga.put(g_a, a)
    ga.sync()
    ijv = np.zeros((m,2), dtype=np.int64)
    random.seed(ga.nodeid()*51 + 1) # different seed for each proc
    for j in range(10):
        itmp = None
        if MIRROR:
            itmp = random.randint(0,lprocs-1)
        else:
            itmp = random.randint(0,nproc-1)
        if itmp == me:
            for loop in range(m):
                ijv[loop,:] = (random.randint(0,n-1),random.randint(0,n-1))
                #if ijv[loop,0] > ijv[loop,1]:
                #    ijv[loop,:] = ijv[loop,::-1] # reverse
            result = ga.gather(g_a, ijv)
            for loop in range(m):
                value = ga.get(g_a, ijv[loop], ijv[loop]+1).flatten()
                if not result[loop] == value:
                    ga.error('gather failed')
    if 0 == me:
        print 'OK'
    ga.destroy(g_a)

def check_scatter(gatype):
    nptype = ga.dtype(gatype)
    if 0 == me:
        print '> Checking scatter (might be slow)...',
    g_a = create_global_array(gatype)
    a = create_local_a(gatype)
    if 0 == me:
        ga.put(g_a, a)
    ga.sync()
    ijv = np.zeros((m,2), dtype=np.int64)
    v = np.zeros(m, dtype=nptype)
    random.seed(ga.nodeid()*51 + 1) # different seed for each proc
    for j in range(10):
        check = None
        if MIRROR:
            check = random.randint(0,lprocs-1) == iproc
        else:
            check = random.randint(0,nproc-1) == me
        if check:
            for loop in range(m):
                ijv[loop,:] = (random.randint(0,n-1),random.randint(0,n-1))
                v[loop] = ijv[loop,0]+ijv[loop,1]
            ga.scatter(g_a, v, ijv)
            for loop in range(m):
                value = ga.get(g_a, ijv[loop], ijv[loop]+1).flatten()
                if not v[loop] == value:
                    ga.error('scatter failed')
    if 0 == me:
        print 'OK'
    ga.destroy(g_a)

def check_print_patch(gatype):
    g_a = create_global_array(gatype)
    a = create_local_a(gatype)
    if n > 7:
        if 0 == me:
            print '> Checking ga.print_patch --- should match'
            print a[2:5,2:7]
        ga.print_patch(g_a, (2,2), (5,7))
    ga.destroy(g_a)

def check_read_inc(gatype):
    if 0 == me:
        print '> Checking ga.read_inc ...',
        print "CHECK NOT IMPLEMENTED" 
    #print 'OK'

def check_fence_and_lock(gatype):
    if 0 == me:
        print '> Checking ga.fence and ga.lock',
    g_a = create_global_array(gatype)
    ga.zero(g_a)
    if not ga.create_mutexes(1):
        ga.error('ga.create_mutexes failed')
    if n < 2:
        ga.error('insufficient n to test ga.fence', n)
    ga.lock(0)
    a = ga.get(g_a) # get original values
    a[:,0] += 1 # add my contribution
    # need to use fence to assure that coms complete before leaving
    # critical section
    ga.init_fence()
    ga.put(g_a, a)
    ga.fence()
    ga.unlock(0)
    if not ga.destroy_mutexes():
        ga.error('mutex not destroyed')
    ga.sync()
    if 0 == me:
        a = ga.get(g_a)
        if not np.all(a[:,0] == nproc):
            ga.error('fence failed')
    if 0 == me:
        print 'OK'

def check(gatype):
    check_zero(gatype)
    check_put_disjoint(gatype)
    check_get(gatype)
    check_accumulate_disjoint(gatype)
    check_accumulate_overlap(gatype)
    check_add(gatype)
    check_dot(gatype)
    check_scale(gatype)
    check_copy(gatype)
    check_gather(gatype)
    check_scatter(gatype)

def check_float():
    check(ga.C_FLOAT)

def check_double():
    check(ga.C_DBL)

def check_complex_float():
    check(ga.C_SCPL)

def check_complex_double():
    check(ga.C_DCPL)

def check_int():
    gatype = ga.C_INT
    check_zero(gatype)
    check_put_disjoint(gatype)
    check_get(gatype)
    check_accumulate_disjoint(gatype)
    check_accumulate_overlap(gatype)
    #check_add(gatype)
    #check_dot(gatype)
    #check_scale(gatype)
    check_copy(gatype)
    check_gather(gatype)
    check_scatter(gatype)
    check_print_patch(gatype)
    check_fence_and_lock(gatype)

def check_long():
    gatype = ga.C_LONG
    check_zero(gatype)
    check_put_disjoint(gatype)
    check_get(gatype)
    check_accumulate_disjoint(gatype)
    check_accumulate_overlap(gatype)
    #check_add(gatype)
    #check_dot(gatype)
    #check_scale(gatype)
    check_copy(gatype)
    check_gather(gatype)
    check_scatter(gatype)
    check_print_patch(gatype)
    check_fence_and_lock(gatype)

def check_gop(nptype):
    if 0 == me:
        print '> checking ga.gop (%s)' % nptype,
    input = np.arange(n, dtype=nptype) + me
    sum = np.arange(n, dtype=nptype)*nproc + (nproc-1)*nproc/2
    output = ga.gop(input, '+')
    if not np.all(output == sum):
        ga.error('ga.gop (%s) error' % nptype)
    if 0 == me:
        print 'OK'

def check_broadcast():
    if 0 == me:
        print '> Checking ga.brdcst',
    buf = [0,0]
    if nproc-1 == me:
        buf = [me,nproc]
    buf = ga.brdcst(buf,nproc-1)
    if buf[0] != nproc-1:
        ga.error('ga.brdcst buf[0] failed')
    if buf[1] != nproc:
        ga.error('ga.brdcst buf[1] failed')
    if 0 == me:
        print 'OK'

def check_wrappers():
    check_gop(np.int32)
    check_gop(np.int64)
    check_gop(np.float32)
    check_gop(np.float64)
    check_gop(np.complex64)
    check_gop(np.complex128)
    check_broadcast()

if __name__ == '__main__':
    main()
