import math
import mathutil
from array import *
from mathutil import printvector, printmatrix, dot
from longfloat import longfloat
from quadrature import Quadrature, QuadratureTest

# Compute the two-scale coefficients for Aplert's multi-wavelets
# This has been separated from the other code since we need extended
# precision to do this.
#
# Import this file using
#
# from twoscalecoeffs import twoscalecoeffs
#
# otherwise you'll end up with a bunch of multiple precision stuff that
# you really don't want

# Redefine a vector to be a list rather than an array of doubles
# so that we can use extended precision

def vector(n):
    return range(n)

def zerovector(n):
    v = vector(n)
    for i in range(n):
        v[i] = longfloat(0)
    return v

def zeromatrix(n,m): 
    a = range(n) 
    for i in range(n): 
        a[i] = zerovector(m) 
    return a 

def pn(x,order):
    '''
    Evaluate the Legendre polynomials up to the given order at x
    defined on [-1,1].
    '''
    p = range(order+1)
    p[0] = 1
    if order == 0:
        return p
    p[1] = x
    for n in range(1,order):
        p[n+1] = n*(x*p[n] - p[n-1])/(n+1) + x*p[n]
    return p

phi_norms = []                          # Cache sqrt(2*n+1)
def phi(x,k):
    '''
    Evaluate the shifted normalized Legendre polynomials up to the
    given order at x defined on [0,1].

    These are also our scaling functions, phi_i(x) , i=0..k-1

    In addition to forming an orthonormal basis on [0,1] we have
    phi_j(1/2-x) = (-1)^j phi_j(1/2+x)

    (the wavelets are similar with phase (-1)^(j+k)).
    '''
    global phi_norms

    if len(phi_norms) == 0:
        for n in range(max(k,100)):
            phi_norms.append(math.sqrt(2*n+1))
    order = k-1
    p = pn(2.*x-1,order)
    for n in range(k):
        p[n] = p[n]*phi_norms[n]  # sqrt(2*n+1)
    return p

def psi(x,k,g0,g1):
    '''
    For debug and demo purposes evaluate the multi-wavelets at x in [0,1]
    given the two-scale coefficients.
    '''
    if x < 0.5:
        p = phi(2.0*x,k)
        g = g0
    else:
        p = phi(2.0*x-1.0,k)
        g = g1

    root2 = math.sqrt(2.0)
    value = mathutil.zerovector(k)
    for i in range(k):
        for j in range(k):
            value[i] = value[i] + root2*g[i][j]*p[j]

    return value


def pnlong(x,order):
    '''
    USE LONGFLOATS
    
    Evaluate the Legendre polynomials up to the given order at x
    defined on [-1,1].
    '''
    p = range(order+1)
    p[0] = longfloat(1)
    x = longfloat(x)
    if order == 0:
        return p
    p[1] = x
    for n in range(1,order):
        p[n+1] = n*(x*p[n] - p[n-1])/(n+1) + x*p[n]
    return p

phi_norms_long = []                          # Cache sqrt(2*n+1)
def philong(x,k):
    '''
    USE LONGFLOATS
    
    Evaluate the shifted normalized Legendre polynomials up to the
    given order at x defined on [0,1].

    These are also our scaling functions, phi_i(x) , i=0..k-1

    In addition to forming an orthonormal basis on [0,1] we have
    phi_j(1/2-x) = (-1)^j phi_j(1/2+x)

    (the wavelets are similar with phase (-1)^(j+k)).
    '''
    global phi_norms_long

    if len(phi_norms_long) == 0:
        for n in range(max(k+1,61)):
            phi_norms_long.append(longfloat(2*n+1).sqrt())
    if k > 60:
        raise "phi_norms_long too small"
    
    order = k-1
    p = pnlong(2*x-1,order)
    for n in range(k):
        p[n] = p[n]*phi_norms_long[n]  # sqrt(2*n+1)
    return p

def h0h1(k):
    '''
    USE LONGFLOATS
    
    Compute the h0 & h1 two-scale coefficients

    h0_ij = sqrt(2)*int(phi_i(x)*phi_j(2x),x=0..1/2)
    h1_ij = sqrt(2)*int(phi_i(x)*phi_j(2x-1),x=1/2..1)

    but
    
    h1_ij = (-1)^(i+j) h0_ij from symmetry property of the phi

    The h0/h1 are useful because

    phi_i(x) = sqrt(2)*sum(j) h0_ij*phi_j(2x) + h1_ij*phi_j(2x-1)

    and
    
    phi_i(2x)   = (1/sqrt(2))*sum(j) h0_ji*phi_j(x) + g0_ji*psi(x)
    phi_i(2x-1) = (1/sqrt(2))*sum(j) h1_ji*phi_j(x) + g1_ji*psi(x)
    '''
    q = Quadrature(k,[longfloat(0),longfloat((-1,1))],uselongfloat=1)
    x = q.points()
    w = q.weights()
    npt = len(x)
    twox = range(npt)
    for i in range(npt):
        twox[i] = x[i] + x[i]

    h0 = zeromatrix(k,k)
    h1 = zeromatrix(k,k)
    
    root2 = longfloat(2).sqrt()
    for t in range(npt):
        p = philong(x[t],k)
        p2= philong(twox[t],k)
        for i in range(k):
            faci = p[i]*w[t]*root2
            for j in range(k):
                h0[i][j] = h0[i][j] + faci*p2[j]

    for i in range(k):
        for j in range(k):
            h1[i][j] = ((-1)**(i+j))*h0[i][j]

    return (h0, h1)

def orthonormalize_wavelets(g0,g1,k):
    '''
    USE LONGFLOATS
    
    Orthonormalize in reverse order
    '''
    for i in range(k-1,-1,-1):
        for attempt in range(1):   # was 2
            # Orthog to previous functions
            for m in range(k-1,i,-1):
                dotim = dot(g0[i]+g1[i],g0[m]+g1[m])
                for j in range(k):
                    g0[i][j] = g0[i][j] - dotim*g0[m][j]
                    g1[i][j] = g1[i][j] - dotim*g1[m][j]
            norm = dot(g0[i]+g1[i],g0[i]+g1[i])
            rnorm = norm.rsqrt()
            for j in range(k):
                g0[i][j] = g0[i][j]*rnorm
                g1[i][j] = g1[i][j]*rnorm

def half_moments(maxm,k):
    '''
    USE LONGFLOATS
    
    Compute the moments of the scaling functions computed
    over half over the interval.

    hmom1[p][j] = int(x^p * phi_j(2*x), x=0..1/2)
    hmom2[p][j] = int((x+0.5)^p * phi_j(2*x), x=0..1/2)

    hmom2 is made from hmom1 using the binomial expansion
    (x+a)^p = a^0*x^p + p*a^1*x^p-1 + p*(p-1)*a^2*x^p-2 / 2 + ... + a^p*x^0

    From Alpert, Beylkin, Gines and Vozovoi

    int(x^p phi_j(x),x=0..1) = sqrt(2*j+1)*(p!)^2 / ((p-j)!(p+j+1)!) p>=j

    Get an extra factor of (1/2)^(m+1) from the change of scale.

    For completeness also record here that
    
    int((x-1)^p phi_j(x),x=0..1) = (-1)^(j+p) *
    .                          sqrt(2*j+1)*(p!)^2 / ((p-j)!(p+j+1)!)

    Moments with p < j are zero due to orthogonality.
    '''

    hmom1 = zeromatrix(maxm+1,k)
    hmom2 = zeromatrix(maxm+1,k)

    factorial = range(maxm+1+k)
    pfac = longfloat(1)
    for p in range(maxm+1+k):
        factorial[p] = pfac
        pfac = pfac * (p+1)

    half = longfloat((-1,1))
    for p in range(maxm+1):
        pfac2 = factorial[p]*factorial[p]
        scale = longfloat(half)**(p+1)
        for j in range(k):
            if p < j:
                hmom1[p][j] = longfloat(0)
            else:
                hmom1[p][j] = scale*(longfloat(2*j+1).sqrt())*pfac2 / \
                             (factorial[p-j]*factorial[p+j+1])

    for p in range(maxm+1):
        for i in range(k):
            fac = longfloat(1)
            for q in range(p+1):
                hmom2[p][i] = hmom2[p][i] + fac*(half**(p-q))*hmom1[q][i]
                fac = fac * (p-q) / (q+1)
    return (hmom1, hmom2)
    

def g0g1(k):
    '''
    USE LONGFLOATS
    
    Compute the g0 & g1 two-scale coefficients which define the
    multi-wavelet in terms of the scaling functions at the finer
    level.

    psi_i(x) = sqrt(2)*sum(j) g0_ij*phi_j(2x) + g1_ij*phi_j(2x-1)

    Projection with sqrt(2)*phi_j(2x) and sqrt(2)*phi_j(2x-1), which
    are the normalized scaling functions on the finer level, yields 

    g0_ij = sqrt(2)*int(psi_i(x)*phi_j(2x),x=0..1/2)
    g1_ij = sqrt(2)*int(psi_i(x)*phi_j(2x-1),x=1/2..1)


    With the original Alpert conditions we have
    
    g1_ij = (-1)^(i+j+k) g0_ij from symmetry properties of the psi

    but you only get this if you impose the condition of additional
    vanishing moments.

    The projection of the additional moments is very numerically
    unstable and forces the use of higher precision for k >= 6.
    '''

    g0 = zeromatrix(k,k)
    g1 = zeromatrix(k,k)

    # Initialize the basis
    rroot2 = longfloat(2).rsqrt()
    for i in range(0,k):
        g0[i][i] =  rroot2
        g1[i][i] = -rroot2

    (h0, h1) = h0h1(k)

    # Project out the moments == orthog to the scaling functions.
    # Need to use the two-scale relation for the scaling functions.
    for i in range(k):
        for attempt in range(1): # was 2
            for m in range(k):
                dotim = dot(g0[i],h0[m])+dot(g1[i],h1[m])
                for j in range(k):
                    g0[i][j] = g0[i][j] - dotim*h0[m][j]
                    g1[i][j] = g1[i][j] - dotim*h1[m][j]

    # Orthog to x^k, x^k+1, ... , x^2k-2
    first = 0
    root2 = longfloat(2).sqrt()
    (hmom1,hmom2) = half_moments(2*k-2,k)

    test = longfloat(1).eps().sqrt()
    for power in range(k,2*k-1):
        # Transform to the wavelet basis
        overlap = zerovector(k)
        for i in range(k):
            for j in range(k):
                overlap[i] = overlap[i] + root2*\
                  (g0[i][j]*hmom1[power][j] + g1[i][j]*hmom2[power][j])

        # Put first function with non-zero overlap in front of the list
        for i in range(first,k):
            if abs(overlap[i]) > test:       # Is there a better test?
                (g0[first],g0[i]) = (g0[i],g0[first])
                (g1[first],g1[i]) = (g1[i],g1[first])
                (overlap[first],overlap[i]) = (overlap[i],overlap[first])
                break

        # Project out the first component as necessary.
        if abs(overlap[first]) > test:       # Is there a better test?
            for i in range(first+1,k):
                scale = overlap[i]/overlap[first]
                for j in range(k):
                    g0[i][j] = g0[i][j] - scale*g0[first][j]
                    g1[i][j] = g1[i][j] - scale*g1[first][j]

        first = first + 1
        
    orthonormalize_wavelets(g0,g1,k)
        
    return (g0, g1)

def fromfile(k):
    '''
    See if the two-scale coeffs have already been computed and stored on disk.
    Return the data are successfully read from disk, otherwise None.
    '''
    fname = '/tmp/two-scale-coeffs.' + str(k)
    try:
        f = open(fname,'rb')
        data = array('d',[])
        data.fromfile(f,4*k*k)
        h0 = mathutil.zeromatrix(k,k)
        h1 = mathutil.zeromatrix(k,k)
        g0 = mathutil.zeromatrix(k,k)
        g1 = mathutil.zeromatrix(k,k)
        pt = 0
        for i in range(k):
            for j in range(k):
                h0[i][j] = data[pt]
                h1[i][j] = data[pt+1]
                g0[i][j] = data[pt+2]
                g1[i][j] = data[pt+3]
                pt = pt + 4
        #print ' Read two-scale coefficients from file',fname
        return (h0,h1,g0,g1)
    except IOError, EOFError:
        return None

def tofile(h0,h1,g0,g1):
    '''
    Write the two-scale coeffs to disk.
    '''
    k = len(h0)
    fname = '/tmp/two-scale-coeffs.' + str(k)
    f = open(fname,'w+b')
    data = array('d',[0]*(4*k*k))
    pt = 0
    for i in range(k):
        for j in range(k):
            data[pt  ] = h0[i][j]
            data[pt+1] = h1[i][j]
            data[pt+2] = g0[i][j]
            data[pt+3] = g1[i][j]
            pt = pt + 4
    data.tofile(f)
    f.close()
    #print ' Wrote two-scale coefficients to file',fname


def twoscalecoeffs(k,usecached=1):
    '''
    Return the two-scale coefficients for Alperts basis of given order

    coarser <= finer
    phi_i(x) = sqrt(2)*sum(j) h0_ij*phi_j(2x) + h1_ij*phi_j(2x-1)
    ditto for psi_i(x) with g instead of h

    s_n,i,l = sum(j) h0_ij*s_n+1,j,2l + h1_ij*s_n+1,j,2l+1
    ditto for d with g instead of h

    finer <= coarser
    phi_i(2x)   = (1/sqrt(2))*sum(j) h0_ji*phi_j(x) + g0_ji*psi(x)
    phi_i(2x-1) = (1/sqrt(2))*sum(j) h1_ji*phi_j(x) + g1_ji*psi(x)

    s_n+1,i,2l   = sum(j) h0_ji*sn,j,l + g0_ji*d_n,j,l
    s_n+1,i,2l+1 = sum(j) h1_ji*sn,j,l + g1_ji*d_n,j,l

    '''
    global phi_norms_long

    if usecached:
        cached = fromfile(k)
        if cached:
            return cached

    nbitssave = longfloat.nbits
    if k < 11:
        longfloat.nbits = 156
    elif k < 19:
        longfloat.nbits = 208
    elif k < 24:
        longfloat.nbits = 260
    elif k < 29:
        longfloat.nbits = 320
    elif k < 33:
        longfloat.nbits = 360
    elif k < 37:
        longfloat.nbits = 400
    elif k < 42:
        longfloat.nbits = 440
    elif k < 46:
        longfloat.nbits = 480
    elif k < 49:
        longfloat.nbits = 520
    else:
        longfloat.nbits = 800

    longfloat.nbits = 800
    print ' Generating two-scale coeffcients for k=%d using %d-bit floating point numbers' \
          % (k, longfloat.nbits)
    
    # Force recomputation of the cached sqrt(2*j+1) values 
    phi_norms_long = []
    
    # Compute the coeffcicients to high precision using longfloats
    (h0l,h1l) = h0h1(k)
    (g0l,g1l) = g0g1(k)

    # Copy them into arrays of floats
    h0 = mathutil.zeromatrix(k,k)
    g0 = mathutil.zeromatrix(k,k)
    h1 = mathutil.zeromatrix(k,k)
    g1 = mathutil.zeromatrix(k,k)

    gerr = 0.0
    herr = 0.0
    for i in range(k):
        for j in range(k):
            hphase = (-1)**(i+j)
            gphase = (-1)**(i+j+k)
            h0[i][j] = float(h0l[i][j])
            h1[i][j] = h0[i][j]*hphase
            g0[i][j] = float(g0l[i][j])
            g1[i][j] = g0[i][j]*gphase
            if abs(h0[i][j]) < 4e-16: h0[i][j] = 0.0
            if abs(h1[i][j]) < 4e-16: h1[i][j] = 0.0
            if abs(g0[i][j]) < 4e-16: g0[i][j] = 0.0
            if abs(g1[i][j]) < 4e-16: g1[i][j] = 0.0
            testerr = float(abs(g0l[i][j]-gphase*g1l[i][j]))
            if g0[i][j] != 0.0: testerr = testerr/abs(g0[i][j])
            gerr = max(gerr,testerr)
            testerr = float(abs(h0l[i][j]-hphase*h1l[i][j]))
            if h0[i][j] != 0.0: testerr = testerr/abs(h0[i][j])
            herr = max(herr,testerr)
            
    print ' maximum symmetry error in h ', float(herr)
    print ' maximum symmetry error in g ', float(gerr)

    longfloat.nbits = nbitssave
    
    if max(float(herr),float(gerr)) > 4e-16:
        raise "twoscalecoeffs have inadequate precision"

    tofile(h0,h1,g0,g1)

    return (h0,h1,g0,g1)

if __name__ == "__main__":

    (h0,h1,g0,g1) = twoscalecoeffs(4)
    print h0
    print h1
    print g0
    print g1
    atop


    print '\n\n Beginning new run \n'

    #QuadratureTest(uselongfloat=1)

    for k in range(1,31):
        print '\n\n Testing multiwavelets k =', k, '\n'

        (h0,h1,g0,g1) = twoscalecoeffs(k)
        # Demonstrate that the scaling functions are ortho-normal
        i = 0
        j = 0
        def f(x):
            global k, i, j
            p = phi(x,k)
            return p[i]*p[j]
        q = Quadrature(k,[0,1.0])
        err = 0.0
        for i in range(k):
            for j in range(i+1):
                overlap = q.integrate(f)
                if i == j: overlap = overlap - 1
                err = max(err,abs(overlap))
        print ' Maximum errror in the overlap of scaling functions', err

        # Demonstrate that the wavelets are ortho-normal
        i = 0
        j = 0
        def f(x):
            global k, i, j, g0, g1
            p = psi(x,k,g0,g1)
            return p[i]*p[j]

        q1 = Quadrature(k,[0,0.5])
        q2 = Quadrature(k,[0.5,1.0])
        err = 0.0
        for i in range(k):
            for j in range(i+1):
                overlap = q1.integrate(f) + q2.integrate(f)
                if i == j: overlap = overlap - 1
                err = max(err,abs(overlap))
        print ' Maximum errror in the overlap of wavelets', err

        # Demonstrate that the wavelets are normal to scaling functions
        i = 0
        j = 0
        def f(x):
            global k, i, j, g0, g1
            p = psi(x,k,g0,g1)
            q = phi(x,k)
            return p[i]*q[j]

        q1 = Quadrature(k,[0,0.5])
        q2 = Quadrature(k,[0.5,1.0])
        err = 0.0
        for i in range(k):
            for j in range(k):
                overlap = q1.integrate(f) + q2.integrate(f)
                err = max(err,abs(overlap))
        print ' Maximum errror in the overlap of wavelets with scaling functions', err

        if k <= 4:
            print '     h0 '
            printmatrix(h0)
            print '     g0 '
            printmatrix(g0)
        
