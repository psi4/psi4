import time
#import numpy as np
import ga.gain as np
from ga.gain import me
import random
import sys

def MC_int(c, s):
    x = np.empty([s], dtype=np.double)
    y = np.empty([s], dtype=np.double)
    sum=0.0
    #np.timer_reset()
    #np.evalflush()
    start=time.time()
    for i in range(c):
        #np.ufunc_random(x,x)
        #np.ufunc_random(y,y)
        x = np.random.random_sample([s])
        y = np.random.random_sample([s])
        np.square(x,x)
        np.square(y,y)
        np.add(x,y,x)
        z = np.less_equal(x, 1)
        sum += np.add.reduce(z)*4.0/s
    sum = sum / c
    #np.evalflush()
    stop=time.time()
    print 'Pi: ', sum, ' with ', s,' samples in sec: ', stop-start,
    print "(Dist) notes: %s"%sys.argv[3]
    #if dist:
    #    print "(Dist) notes: %s"%sys.argv[4]
    #else:
    #    print "(Non-Dist) notes: %s"%sys.argv[4]


def main():
    #D=int(sys.argv[1])
    N=int(sys.argv[1])
    C=int(sys.argv[2])
    MC_int(C, N)

if __name__ == '__main__':
    #import cProfile
    #cProfile.run("main()", "profile.%s.prof" % str(me()))
    main()
