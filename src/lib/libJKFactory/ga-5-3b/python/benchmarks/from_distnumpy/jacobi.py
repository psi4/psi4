import time
import sys
#import numpy as np
import ga.gain as np
from ga.gain import me

def jacobi(A, B, tol=0.005, forcedIter=0):
    '''itteratively solving for matrix A with solution vector B
       tol = tolerance for dh/h
       init_val = array of initial values to use in the solver
    '''
    h = np.zeros(np.shape(B), float)
    dmax = 1.0
    n = 0
    tmp0 = np.empty(np.shape(A), float)
    tmp1 = np.empty(np.shape(B), float)
    AD = np.diagonal(A)
    #np.timer_reset()
    #np.evalflush()
    t1 = time.time()
    while (forcedIter and forcedIter > n) or \
          (forcedIter == 0 and dmax > tol):
        n += 1
        np.multiply(A,h,tmp0)
        np.add.reduce(tmp0,1,out=tmp1)
        tmp2 = AD
        np.subtract(B, tmp1, tmp1)
        np.divide(tmp1, tmp2, tmp1)
        hnew = h + tmp1
        np.subtract(hnew,h,tmp2)
        np.divide(tmp2,h,tmp1)
        np.absolute(tmp1,tmp1)
        dmax = np.maximum.reduce(tmp1)
        h = hnew
        print dmax


    #np.evalflush()
    t1 = time.time() - t1

    print 'Iter: ', n, ' size: ', np.shape(B),' time: ', t1,
    print "(Dist) notes: %s"%sys.argv[4]
    #if A.dist():
    #    print "(Dist) notes: %s"%sys.argv[4]
    #else:
    #    print "(Non-Dist) notes: %s"%sys.argv[4]


    return h

def main():
    d = int(sys.argv[1])
    size = int(sys.argv[2])
    iter = int(sys.argv[3])

    #A = array([[4, -1, -1, 0], [-1, 4, 0, -1], [-1, 0, 4, -1], [0, -1, -1, 4]], float, dist=d)
    #B = array([1,2,0,1], float, dist=d)

    #A = np.zeros([size,size], dtype=float)
    #np.ufunc_random(A,A)
    A = np.random.random_sample([size,size])

    #B = np.zeros([size], dtype=float)
    #np.ufunc_random(B,B)
    B = np.random.random_sample([size])

    C = jacobi(A, B, 0.10, forcedIter=iter)

if __name__ == '__main__':
    #import cProfile
    #cProfile.run("main()", "profile.%s.prof" % str(me()))
    main()
