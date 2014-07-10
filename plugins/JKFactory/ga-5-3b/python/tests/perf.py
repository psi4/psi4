"""
perf.py is used to test performance of put, get and accumulate.
It has to be executed on four processors.
Remote operations access data on processors 1,2,3 in the round-robin way.
"""
import time
import sys

import mpi4py.MPI
from ga4py import ga

import numpy as np

nproc = ga.nnodes()
me = ga.nodeid()

def main():
    if 4 != nproc and 0 == me:
        ga.error('Program requires 4 GA processes; nproc=%s' % nproc)
    test2D()
    test1D()
    if 0 == me:
        print 'All tests successful'

def test2D():
    n = 1024
    buf = np.zeros((n,n), dtype=np.float64)
    chunk = np.asarray([1,3,4,9,16,24,30,48,64,91,128,171,256,353,440,512])
    g_a = ga.create(ga.C_DBL, (n,n), 'a')
    if 0 == g_a:
        ga.error('ga.create failed')
    buf[:] = 0.01
    ga.zero(g_a)
    if 0 == me:
        print (' Performance of GA get, put & acc'
                ' for square sections of array[%d,%d]' % (n,n))
    lo,hi = ga.distribution(g_a, me)
    # local ops
    TestPutGetAcc(g_a, n, chunk, buf, lo, hi, True)
    # remote ops
    TestPutGetAcc(g_a, n, chunk, buf, lo, hi, False)

def TestPutGetAcc(g_a, n, chunk, buf, lo, hi, local):
    if 0 == me:
        print ''
        if local:
            print 'Local 2-D Array Section'
        else:
            print 'Remote 2-D Array Section'
        print '%15s %19s %19s %19s' % ('section', 'get', 'put', 'accumulate')
        print '%7s %7s %9s %9s %9s %9s %9s %9s' % (
                'bytes','dim','usec','MB/s','usec','MB/s','usec','MB/s')
    ga.sync()
    bytes = 0
    jump = 0
    num_chunks = len(chunk)
    for loop in range(num_chunks):
        tg = 0.0
        tp = 0.0
        ta = 0.0
        bytes = 8*chunk[loop]*chunk[loop] # how much data is accessed
        jump = n/(60*(loop+1)) # jump between consecutive patches
        if loop+1 == num_chunks:
            jump = 0
        # everybody touches own data
        ga.fill(g_a, me*(loop+1), [0,0], [n,n])
        if 0 == me:
            tg = time_get(g_a, lo, hi, buf, chunk[loop], jump, local)
        else:
            time.sleep(1)
        # everybody touches own data
        ga.fill(g_a, me*(loop+1), [0,0], [n,n])
        if 0 == me:
            tp = time_put(g_a, lo, hi, buf, chunk[loop], jump, local)
        else:
            time.sleep(1)
        # everybody touches own data
        ga.fill(g_a, me*(loop+1), [0,0], [n,n])
        if 0 == me:
            ta = time_acc(g_a, lo, hi, buf, chunk[loop], jump, local)
        else:
            time.sleep(1)
        if 0 == me:
            print '%7d %7d %8.3e %8.3e %8.3e %8.3e %8.3e %8.3e' % (
                    bytes, chunk[loop],
                    tg/1e-6, 1e-6*bytes/tg,
                    tp/1e-6, 1e-6*bytes/tp,
                    ta/1e-6, 1e-6*bytes/ta)

def time_get(g_a, lo, hi, buf, chunk, jump, local):
    count = 0
    rows = hi[0]-lo[0]
    cols = hi[1]-lo[1]
    shifti = [rows, 0, rows]
    shiftj = [0, cols, cols]
    seconds = time.time()
    # distance between consecutive patches increased by jump
    # to destroy locality of reference
    for ilo in range(lo[0], hi[0]-chunk-jump+1, chunk+jump):
        ihi = ilo + chunk
        for jlo in range(lo[1], hi[1]-chunk-jump+1, chunk+jump):
            jhi = jlo + chunk
            count += 1
            if local:
                llo = [ilo,jlo]
                lhi = [ihi,jhi]
                ignore = ga.get(g_a, llo, lhi, buf[ga.zip(llo,lhi)])
            else:
                index = count%3
                llo = [ilo+shifti[index],jlo+shiftj[index]]
                lhi = [ihi+shifti[index],jhi+shiftj[index]]
                ignore = ga.get(g_a, llo, lhi, buf[ilo:ihi,jlo:jhi])
    seconds = time.time() - seconds
    return seconds/count

def time_put(g_a, lo, hi, buf, chunk, jump, local):
    count = 0
    rows = hi[0]-lo[0]
    cols = hi[1]-lo[1]
    shifti = [rows, 0, rows]
    shiftj = [0, cols, cols]
    seconds = time.time()
    # distance between consecutive patches increased by jump
    # to destroy locality of reference
    for ilo in range(lo[0], hi[0]-chunk-jump+1, chunk+jump):
        ihi = ilo + chunk
        for jlo in range(lo[1], hi[1]-chunk-jump+1, chunk+jump):
            jhi = jlo + chunk
            count += 1
            if local:
                llo = [ilo,jlo]
                lhi = [ihi,jhi]
                ga.put(g_a, buf[ga.zip(llo,lhi)], llo, lhi)
            else:
                index = count%3
                llo = [ilo+shifti[index],jlo+shiftj[index]]
                lhi = [ihi+shifti[index],jhi+shiftj[index]]
                ga.put(g_a, buf[ilo:ihi,jlo:jhi], llo, lhi)
    seconds = time.time() - seconds
    return seconds/count

def time_acc(g_a, lo, hi, buf, chunk, jump, local):
    count = 0
    rows = hi[0]-lo[0]
    cols = hi[1]-lo[1]
    shifti = [rows, 0, rows]
    shiftj = [0, cols, cols]
    seconds = time.time()
    # distance between consecutive patches increased by jump
    # to destroy locality of reference
    for ilo in range(lo[0], hi[0]-chunk-jump+1, chunk+jump):
        ihi = ilo + chunk
        for jlo in range(lo[1], hi[1]-chunk-jump+1, chunk+jump):
            jhi = jlo + chunk
            count += 1
            if local:
                llo = [ilo,jlo]
                lhi = [ihi,jhi]
                ga.acc(g_a, buf[ga.zip(llo,lhi)], llo, lhi, 1)
            else:
                index = count%3
                llo = [ilo+shifti[index],jlo+shiftj[index]]
                lhi = [ihi+shifti[index],jhi+shiftj[index]]
                ga.acc(g_a, buf[ilo:ihi,jlo:jhi], llo, lhi, 1)
    seconds = time.time() - seconds
    return seconds/count

def test1D():
    n = 1024*1024
    buf = np.zeros(n/4, dtype=np.float64)
    chunk = np.asarray([1,9,16,81,256,576,900,2304,4096,8281,
        16384,29241,65536,124609,193600,262144])
    g_a = ga.create(ga.C_DBL, (n,), 'a')
    if 0 == g_a:
        ga.error('ga.create failed')
    buf[:] = 0.01
    ga.zero(g_a)
    if 0 == me:
        print ''
        print ''
        print ''
        print (' Performance of GA get, put & acc'
                ' for 1-dimensional sections of array[%d]' % n)
    lo,hi = ga.distribution(g_a, me)
    # local ops
    TestPutGetAcc1(g_a, n, chunk, buf, lo, hi, True)
    # remote ops
    TestPutGetAcc1(g_a, n, chunk, buf, lo, hi, False)

def TestPutGetAcc1(g_a, n, chunk, buf, lo, hi, local):
    if 0 == me:
        print ''
        if local:
            print 'Local 1-D Array Section'
        else:
            print 'Remote 1-D Array Section'
        print '%15s %19s %19s %19s' % ('section', 'get', 'put', 'accumulate')
        print '%7s %7s %9s %9s %9s %9s %9s %9s' % (
                'bytes','dim','usec','MB/s','usec','MB/s','usec','MB/s')
    ga.sync()
    bytes = 0
    jump = 0
    num_chunks = len(chunk)
    for loop in range(num_chunks):
        tg = 0.0
        tp = 0.0
        ta = 0.0
        bytes = 8*chunk[loop] # how much data is accessed
        jump = n/(6000*(loop+1)) # jump between consecutive patches
        if loop+1 == num_chunks:
            jump = 0
        # everybody touches own data
        ga.fill(g_a, me*(loop+1), [0], [n])
        if 0 == me:
            tg = time_get1(g_a, lo, hi, buf, chunk[loop], jump, local)
        else:
            time.sleep(1)
        # everybody touches own data
        ga.fill(g_a, me*(loop+1), [0], [n])
        if 0 == me:
            tp = time_put1(g_a, lo, hi, buf, chunk[loop], jump, local)
        else:
            time.sleep(1)
        # everybody touches own data
        ga.fill(g_a, me*(loop+1), [0], [n])
        if 0 == me:
            ta = time_acc1(g_a, lo, hi, buf, chunk[loop], jump, local)
        else:
            time.sleep(1)
        if 0 == me:
            print '%7d %7d %8.3e %8.3e %8.3e %8.3e %8.3e %8.3e' % (
                    bytes, chunk[loop],
                    tg/1e-6, 1e-6*bytes/tg,
                    tp/1e-6, 1e-6*bytes/tp,
                    ta/1e-6, 1e-6*bytes/ta)

def time_get1(g_a, lo, hi, buf, chunk, jump, local):
    count = 0
    rows = hi[0]-lo[0]
    shift = [3*rows, 2*rows, rows]
    seconds = time.time()
    # distance between consecutive patches increased by jump
    # to destroy locality of reference
    for ilo in range(lo[0], hi[0]-chunk-jump+1, chunk+jump):
        ihi = ilo+chunk
        count += 1
        if local:
            ignore = ga.get(g_a, [ilo], [ihi], buf[ilo:ihi])
        else:
            index = count%3
            llo = ilo+shift[index]
            lhi = ihi+shift[index]
            ignore = ga.get(g_a, llo, lhi, buf[ilo:ihi])
    seconds = time.time() - seconds
    return seconds/count

def time_put1(g_a, lo, hi, buf, chunk, jump, local):
    count = 0
    rows = hi[0]-lo[0]
    shift = [rows, 2*rows, 3*rows]
    seconds = time.time()
    # distance between consecutive patches increased by jump
    # to destroy locality of reference
    for ilo in range(lo[0], hi[0]-chunk-jump+1, chunk+jump):
        ihi = ilo+chunk
        count += 1
        if local:
            ga.put(g_a, buf[ilo:ihi], [ilo], [ihi])
        else:
            index = count%3
            llo = ilo+shift[index]
            lhi = ihi+shift[index]
            ga.put(g_a, buf[ilo:ihi], llo, lhi)
    seconds = time.time() - seconds
    return seconds/count

def time_acc1(g_a, lo, hi, buf, chunk, jump, local):
    # Note: differs from test.F because the passed buffer must be the same
    # size/shape as the patch. The slicing should be fast as the buffer is 1D
    # and contiguous (and so is the slice).
    count = 0
    rows = hi[0]-lo[0]
    shift = [rows, 2*rows, 3*rows]
    seconds = time.time()
    # distance between consecutive patches increased by jump
    # to destroy locality of reference
    for ilo in range(lo[0], hi[0]-chunk-jump+1, chunk+jump):
        ihi = ilo+chunk
        count += 1
        if local:
            ga.acc(g_a, buf[ilo:ihi], [ilo], [ihi], 1.0)
        else:
            index = count%3
            ga.acc(g_a, buf[ilo:ihi], ilo+shift[index], ihi+shift[index], 1.0)
    seconds = time.time() - seconds
    return seconds/count

if __name__ == '__main__':
    main()
