#!/usr/bin/env python

"""
This script compares different ways of implementing an iterative
procedure to solve Laplace's equation.  These provide a general
guideline to using Python for high-performance computing and also
provide a simple means to compare the computational time taken by the
different approaches.  The script compares functions implemented in
pure Python, Numeric, weave.blitz, weave.inline, fortran (via f2py)
and Pyrex.  The function main(), additionally accelerates the pure
Python version using Psyco and provides some numbers on how well that
works.  To compare all the options you need to have Numeric, weave,
f2py, Pyrex and Psyco installed.  If Psyco is not installed the script
will print a warning but will perform all other tests.

The fortran and pyrex modules are compiled using the setup.py script
that is provided with this file.  You can build them like so:

  python setup.py build_ext --inplace


Author: Prabhu Ramachandran <prabhu_r at users dot sf dot net>
License: BSD
Last modified: Sep. 18, 2004
"""

#import numpy
import ga
import ga.gain as numpy
import sys, time

class PrintZero(object):
    def __init__(self):
        self.me = ga.nodeid()
        self.stdout = sys.stdout
    def write(self, something):
        if not self.me:
            self.stdout.write(something)
    def flush(self):
        if not self.me:
            self.stdout.flush()
sys.stdout = PrintZero()

def timer():
    #return time.clock()
    return time.time()

class Grid:
    
    """A simple grid class that stores the details and solution of the
    computational grid."""
    
    def __init__(self, nx=10, ny=10, xmin=0.0, xmax=1.0,
                 ymin=0.0, ymax=1.0):
        self.xmin, self.xmax, self.ymin, self.ymax = xmin, xmax, ymin, ymax
        self.dx = float(xmax-xmin)/(nx-1)
        self.dy = float(ymax-ymin)/(ny-1)
        self.u = numpy.zeros((nx, ny), 'd')
        # used to compute the change in solution in some of the methods.
        self.old_u = self.u.copy()        

    def setBC(self, l, r, b, t):        
        """Sets the boundary condition given the left, right, bottom
        and top values (or arrays)"""        
        self.u[0, :] = l
        self.u[-1, :] = r
        self.u[:, 0] = b
        self.u[:,-1] = t
        self.old_u = self.u.copy()

    def setBCFunc(self, func):
        """Sets the BC given a function of two variables."""
        xmin, ymin = self.xmin, self.ymin
        xmax, ymax = self.xmax, self.ymax
        x = numpy.arange(xmin, xmax + self.dx*0.5, self.dx)
        y = numpy.arange(ymin, ymax + self.dy*0.5, self.dy)
        self.u[0 ,:] = func(xmin,y)
        self.u[-1,:] = func(xmax,y)
        self.u[:, 0] = func(x,ymin)
        self.u[:,-1] = func(x,ymax)

    def computeError(self):        
        """Computes absolute error using an L2 norm for the solution.
        This requires that self.u and self.old_u must be appropriately
        setup."""        
        # use flat
        v = (self.u - self.old_u).flat
        return numpy.sqrt(numpy.dot(v,v))
        # don't use flat
        #v = (self.u - self.old_u)
        #return (v*v).sum()
    

class LaplaceSolver:
    
    """A simple Laplacian solver that can use different schemes to
    solve the problem."""
    
    def __init__(self, grid, stepper='numeric'):
        self.grid = grid
        self.setTimeStepper(stepper)

    def slowTimeStep(self, dt=0.0):
        """Takes a time step using straight forward Python loops."""
        g = self.grid
        nx, ny = g.u.shape        
        dx2, dy2 = g.dx**2, g.dy**2
        dnr_inv = 0.5/(dx2 + dy2)
        # small change here so we can see gain vs numpy slow
        #u = g.u
        u = g.u.get()

        err = 0.0
        for i in range(1, nx-1):
            for j in range(1, ny-1):
                tmp = u[i,j]
                u[i,j] = ((u[i-1, j] + u[i+1, j])*dy2 +
                          (u[i, j-1] + u[i, j+1])*dx2)*dnr_inv
                diff = u[i,j] - tmp
                err += diff*diff

        return numpy.sqrt(err)
        
    def numericTimeStep(self, dt=0.0):
        """Takes a time step using a numeric expressions."""
        g = self.grid
        dx2, dy2 = g.dx**2, g.dy**2
        dnr_inv = 0.5/(dx2 + dy2)
        u = g.u
        g.old_u = u.copy()

        # The actual iteration
        u[1:-1, 1:-1] = ((u[0:-2, 1:-1] + u[2:, 1:-1])*dy2 + 
                         (u[1:-1,0:-2] + u[1:-1, 2:])*dx2)*dnr_inv

        ## first op, assign directly to result
        #numpy.multiply(g.old_u[0:-2, 1:-1], dy2, u[1:-1, 1:-1])
        ## creates a temporary array
        #tmp = g.old_u[2:, 1:-1]*dy2
        ## accumulate temporary array into result
        #numpy.add(tmp, u[1:-1, 1:-1], u[1:-1, 1:-1])
        ## does NOT create a temporary array since we have one already
        #numpy.multiply(g.old_u[1:-1,0:-2], dx2, tmp)
        ## accumulate temporary array into result
        #numpy.add(tmp, u[1:-1, 1:-1], u[1:-1, 1:-1])
        ## does NOT create a temporary array since we have one already
        #numpy.multiply(g.old_u[1:-1, 2:], dx2, tmp)
        ## accumulate temporary array into result
        #numpy.add(tmp, u[1:-1, 1:-1], u[1:-1, 1:-1])
        ## last operation, replace result
        #numpy.multiply(dnr_inv, u[1:-1, 1:-1], u[1:-1, 1:-1])

        return g.computeError()

    def setTimeStepper(self, stepper='numeric'):        
        """Sets the time step scheme to be used while solving given a
        string which should be one of ['slow', 'numeric']."""
        if stepper == 'slow':
            self.timeStep = self.slowTimeStep
        elif stepper == 'numeric':
            self.timeStep = self.numericTimeStep
        else:
            self.timeStep = self.numericTimeStep            
                
    def solve(self, n_iter=0, eps=1.0e-16):        
        """Solves the equation given an error precision -- eps.  If
        n_iter=0 the solving is stopped only on the eps condition.  If
        n_iter is finite then solution stops in that many iterations
        or when the error is less than eps whichever is earlier.
        Returns the error if the loop breaks on the n_iter condition
        and returns the iterations if the loop breaks on the error
        condition."""        
        err = self.timeStep()
        count = 1

        while err > eps:
            if n_iter and count >= n_iter:
                return err
            err = self.timeStep()
            count = count + 1
            print count

        return count


def BC(x, y):    
    """Used to set the boundary condition for the grid of points.
    Change this as you feel fit."""    
    return (x**2 - y**2)


def test(nmin=5, nmax=30, dn=5, eps=1.0e-16, n_iter=0, stepper='numeric'):
    iters = []
    n_grd = numpy.arange(nmin, nmax, dn)
    times = []
    for i in n_grd:
        g = Grid(nx=i, ny=i)
        g.setBCFunc(BC)
        s = LaplaceSolver(g, stepper)
        t1 = timer()
        iters.append(s.solve(n_iter=n_iter, eps=eps))
        dt = timer() - t1
        times.append(dt)
        print "Solution for nx = ny = %d, took %f seconds"%(i, dt)
    return (n_grd**2, iters, times)


def time_test(nx=500, ny=500, eps=1.0e-16, n_iter=100, stepper='numeric'):
    g = Grid(nx, ny)
    g.setBCFunc(BC)
    s = LaplaceSolver(g, stepper)
    t = timer()
    s.solve(n_iter=n_iter, eps=eps)
    return timer() - t
    

def main(n=1000, n_iter=100):
    print "Doing %d iterations on a %dx%d grid"%(n_iter, n, n)
    i = 'numeric'
    print i,
    sys.stdout.flush()
    print "took", time_test(n, n, stepper=i, n_iter=n_iter), "seconds"

    #print "slow (1 iteration)",
    #sys.stdout.flush()
    #s = time_test(n, n, stepper='slow', n_iter=1)
    #print "took", s, "seconds"
    #print "%d iterations should take about %f seconds"%(n_iter, s*n_iter)

if __name__ == "__main__":
    if False:
        import cProfile
        me = ga.nodeid()
        cProfile.run("main()", "profile_%s.prof" % me)
    else:
        main()
    #if me == 0:
    #    ga.print_stats()
