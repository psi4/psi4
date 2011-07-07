
# The stuff in this file is now obsolete and is
# here to keep the old twoscalecoeff code working.

# Don't write new code using this stuff.


from array import * 
from math import sqrt 

def vector(n):
    return array('d',[0.0]*n)

def zerovector(n): 
    return vector(n)

def zeromatrix(n,m): 
    a = range(n) 
    for i in range(n): 
        a[i] = zerovector(m) 
    return a 

def copyvector(x):
    n = len(x)
    a = vector(n)
    for i in range(n): 
        a[i] = x[i] 
    return a 

def copymatrix(x): 
    n = len(x) 
    m = len(x[0]) 
    a = zeromatrix(n,m) 
    for i in range(n): 
        for j in range(m): 
            a[i][j] = x[i][j] 
    return a 

def transpose(x): 
    n = len(x) 
    m = len(x[0]) 
    a = zeromatrix(m,n) 
    for i in range(n): 
        for j in range(m): 
            a[j][i] = x[i][j] 
    return a 

def dot(a,b): 
    sum = 0
    for i in range(len(a)): 
        sum = sum + a[i]*b[i] 
    return sum 

def vector_norm(x):
    '''
    A translation of the BLAS routine dnrm2
    '''
    n = len(x)
    if n < 1:
        return 0.0
    elif n == 1:
        return abs(x[0])
    else:
        scale = 0.0
        ssq = 1.0
        for i in range(n):
            if x[i] != 0.0:
                absxi = abs(x[i])
                if scale < absxi:
                    ssq = 1.0 + ssq*(scale/absxi)**2
                    scale = absxi
                else:
                    ssq = ssq + (absxi/scale)**2
        return scale*sqrt(ssq)

def mxm(a,b): 
    n = len(a) 
    k = len(a[0]) 
    kk= len(b) 
    m = len(b[0]) 
    if kk != k: 
        raise "matrices do not conform for multiplication" 
    c = zeromatrix(n,m) 
    for i in range(n): 
        for j in range(m): 
            sum = 0.0 
            for l in range(k): 
                sum = sum + a[i][l]*b[l][j] 
            c[i][j] = sum 
    return c

sparse_mtxm_nflop = 0.0
def sparse_mtxm(a,b):
    global sparse_mtxm_nflop
    '''
    C = AT * B using sparsity in both.  Currently assumes full storage.
    '''
    k = len(a) 
    kk= len(b) 
    n = len(a[0]) 
    m = len(b[0]) 
    if kk != k: 
        raise "matrices do not conform for multiplication" 
    ct = zeromatrix(m,n)  # This is the transpose of the result
    nflop = 0.0
    for l in range(k):
        al = a[l]
        bl = b[l]
        list = [] # List of non-zero a[l][i] for given l
        for i in range(n):
            if al[i]: list.append(i)
        numl = float(len(list))
        for j in range(m):
            blj = bl[j]
            if blj:
                ctj = ct[j]
                nflop = nflop + numl
                for i in list:
                    ctj[i] = ctj[i] + al[i]*blj
    sparse_mtxm_nflop = sparse_mtxm_nflop + nflop
    return transpose(ct)

def mxv(a,b): 
    n = len(a) 
    k = len(a[0]) 
    kk = len(b) 
    if k != kk: 
        raise "matrix and vector do not conform for multiplication" 
    c = zerovector(n) 
    for i in range(n): 
        c[i] = dot(a[i],b) 
    return c 

def printvector(a): 
    n = len(a)
    for i in range(n): 
        print ("%+15.8e "%float(a[i])),
        if ((i+1)%5) == 0 and i != (n-1):
            print " "
    print " " 

def printmatrix(a): 
    n = len(a)
    m = len(a[0])
    for jlo in range(0,m,5):
        print "    ",
        for j in range(jlo,min(jlo+5,m)):
            print ("     %5i       " % j),
        print " "
        for i in range(n):
            print ("%4i  " % i),
            for j in range(jlo,min(jlo+5,m)):
                print ("%+15.8e " % float(a[i][j])),
            print " "

def numderiv(func,x,step,eps): 
    ''' 
    Use central differences to compute the gradient and diagonal 
    elements of the Hessian. 
    func(x) = function to be differentiated 
    x[] = point at which to differentiate 
    step[] = remembers finite difference step between 
    .        successive calls.  Set to zero on first call 
    .        or set close to appropriate value 
    eps = expected precision in func 
    
    Some care is taken to adjust the step so that the gradient and 
    Hessian diagonal are estimated with about 4 digits of precision 
    but some noise is unavaoidable due either to the noise in the 
    function or cubic/higher terms in the Taylor expansion. 
    ''' 
    
    n = len(x) 
    g = zerovector(n) 
    h = zerovector(n) 
    f0 = func(x) 
    for i in range(n): 
        if step[i] == 0.0: 
            step[i] = max(abs(x[i])*0.01,0.0001) 
        xi = x[i] 
        while 1: 
            x[i] = xi + step[i] 
            f1 = func(x)
            measure = 1e4*eps*max(abs(f0),eps)
            if abs(f1-f0) < measure:
                #print ' Increasing step ',i,step[i],abs(f1-f0) 
                step[i] = step[i]*2.0 
            elif abs(f1-f0) > 10.0*measure:
                #print ' Decreasing step ',i,step[i],abs(f1-f0) 
                step[i] = step[i]/3.0 
            else: 
                break 
        x[i] = xi - step[i] 
        fm1 = func(x) 
        x[i] = xi 
        g[i] = (f1 - fm1)/(2*step[i]) 
        h[i] = (f1 + fm1 - 2.0*f0)/(step[i]*step[i]) 

    return (f0,g,h) 



def quadfit(alpha0, f0, alpha1, f1, alpha2, f2): 
    ''' 
    Given 3 points compute the gradient and hessian at point 0 
    using a quadratic fit. 
    ''' 
    delta1 = alpha1 - alpha0 
    delta2 = alpha2 - alpha0 
    d1 = (f1 - f0)/delta1 
    d2 = (f2 - f0)/delta2 
    h0 = 2.0*(d1 - d2)/(delta1-delta2) 
    g0 = d1 - 0.5*h0*delta1 

    test1 = f0 + g0*delta1 + 0.5*h0*delta1*delta1 
    test2 = f0 + g0*delta2 + 0.5*h0*delta2*delta2 
    return (f0, g0, h0) 

def takestep(x0, s, alpha): 
    x = zerovector(len(x0)) 
    for j in range(len(x)): 
        x[j] = x0[j] + s[j]*alpha 
    return x 

def quadratic_step(trust, g0, h0, prnt=1): 
    if h0 > 0: 
        delta2 = -g0/h0 
        if abs(delta2) > trust:
            if prnt:
                print "                    Step restriction: %f %f " % (delta2, trust) 
            delta2 = abs(trust*delta2)/delta2 
    else:
        if prnt:
            print "                    Negative curvature " 
        delta2 = -abs(trust*g0)/g0 
    return delta2 



def linesearch(func, x0, s, lsgrad, eps, prnt=1): 
    # Assume here that some halfway OK preconditioning 
    # is being used so we expect a step around unity. 
    # Nevertheless, must exercise some caution. 

    # First step in small increments until we've either 
    # bracketed the minimum or gone downhil with enough 
    # energy difference to start fitting

    if prnt:
        print " Line search: step   alpha     grad     hess        value" 
        print "              ---- --------- -------- --------  ----------------" 

    trust = 0.2 
    alpha0 = 0.0 
    f0 = func(x0)
    if prnt:
        print "                   %9.2e %8.1e          %16.8e" % \
              (alpha0, lsgrad, f0) 

    if lsgrad < 0: 
        alpha1 = alpha0 + trust 
    else: 
        alpha1 = alpha0 - trust 
    f1 = func(takestep(x0,s,alpha1))
    if prnt:
        print "                   %9.2e                   %16.8e" % \
              (alpha1, f1) 

    snorm = vector_norm(s)
    
    while f1 > f0: 
        if trust*snorm < sqrt(eps):
            if prnt:
                print " system is too badly conditioned for initial step" 
            return (alpha0,f0)          # Cannot seem to find my way 
        trust = trust * 0.5 
        if lsgrad < 0: 
            alpha1 = alpha0 + trust 
        else: 
            alpha1 = alpha0 - trust 
        f1 = func(takestep(x0,s,alpha1))
        if prnt:
            print "                   %9.2e                   %16.8e" % \
                  (alpha1, f1) 
        
    g0 = lsgrad 
    h0 = (f1-f0-alpha1*g0)/alpha1**2 
    if f1 < f0: 
        g0 = g0 + h0*(alpha1 - alpha0) 
        alpha0, alpha1, f0, f1 = alpha1, alpha0, f1, f0 
    
    alpha2 = alpha0 + quadratic_step(trust,g0,h0,prnt) 

    nbackstep =0 

    for iter in range(1,10): 
        f2 = func(takestep(x0,s,alpha2)) 
        #print ' alphas ', alpha0, alpha1, alpha2 
        #print ' fs     ', f0, f1, f2 
        
        if iter == 1: 
            f2prev = f2 
            
        # Check for convergence or insufficient precision to proceed further 
        if (abs(f0-f1)<(10*eps)) and (abs(f1-f2)<(10*eps)):
            if prnt:
                print "                  ", 
                print " Insufficient precision ... terminating LS" 
            break 
        
        if (f2-f2prev) > 0: 
            # New point is higher than previous worst 
            if  nbackstep < 3: 
                nbackstep = nbackstep + 1
                if prnt:
                    print "                  ", 
                    print " Back stepping due to uphill step" 
                trust = max(0.01,0.2*abs(alpha2 - alpha0))  # Reduce trust radius 
                alpha2 = alpha0 + 0.2*(alpha2 - alpha0) 
                continue 
        elif (f2-f0) < 0: 
            trust = min(4.0,trust*2.0)  # Seem to have narrowed the search 

        nbackstep = 0 

        f2prev = f2 

        # Order in increasing energy 
        if f1 < f0: 
            alpha0, alpha1, f0, f1 = alpha1, alpha0, f1, f0 
        if f2 < f0: 
            alpha0, alpha2, f0, f2 = alpha2, alpha0, f2, f0 
        if f2 < f1: 
            alpha1, alpha2, f1, f2 = alpha2, alpha1, f2, f1 

        (f0, g0, h0) = quadfit(alpha0, f0, alpha1, f1, alpha2, f2)
        if prnt:
            print "              %4i %9.2e %8.1e %8.1e %16.8e" % \
                  (iter, alpha0, g0, h0, f0) 

        if (h0>0.0) and (abs(g0) < 0.03*abs(lsgrad)):
            if prnt:
                print "                  ", 
                print " gradient reduced 30-fold ... terminating LS" 
            break 
        
        # Determine the next step 
        delta = quadratic_step(trust,g0,h0,prnt) 
        alpha2 = alpha0 + delta 
        df = g0*delta + 0.5*h0*delta*delta 
        if abs(df) < 10.0*eps:
            if prnt:
                print "                  ", 
                print " projected energy reduction < 10*eps ... terminating LS" 
            break 

    return (alpha0, f0) 

def jacobi(a): 
    ''' 
    Diagonalize a real symmetric matrix using the variable threshold 
    cyclic Jacobi method.   The input matrix is unmodified.
    
    (v,e) = jacobi(a) 
    
    Input: a[n][m] is a real symmetric matrix 

    Returns: (v,e) where v is the list of eigenvectors and e is an 
    vector of the corresponding eigenvalues in ascending order. 
    v[k] is a vector containing the kth eigenvector.  These satisfy 

    A*Vt = Vt*e 

    or 

    V*A = e*V 

    or 
    
    sum(j)(a[i][j]v[k][j]) = e[k]*v[k][i] 
    ''' 
    a = copymatrix(a) 
    n = len(a) 
    m = len(a[0]) 
    if n != m: 
        raise 'Matrix must be square' 
    for i in range(n): 
        for j in range(m): 
            if a[i][j] != a[j][i]: 
                raise ' Matrix must be symmetric' 

    tolmin = 5.0e-16 
    tol = 1e-2 

    v = zeromatrix(n,n) 
    for i in range(n): 
        v[i][i] = 1.0 

    maxd = 0.0 
    for i in range(n): 
        maxd = max(abs(a[i][i]),maxd) 

    nrotsum = 0
    for iter in range(50):
        maxdaij = 0.0
        nrot = 0 
        for i in range(n): 
            for j in range(i+1,n):  # j>i
                aii = a[i][i] 
                ajj = a[j][j]
                aij = a[i][j]
                daij = abs(aij)
                maxdaij = max(maxdaij,daij/maxd)
                if daij > tol*maxd: # Screen small elements 
                    s = ajj - aii 
                    ds = abs(s) 
                    if daij > (tolmin*ds): # Check for sufficient precision
                        nrot = nrot + 1 
                        if (tolmin*daij) > ds: 
                            c = s = 1/sqrt(2.) 
                        else: 
                            t = aij/s 
                            u = 0.25/sqrt(0.25+t*t) 
                            c = sqrt(0.5+u) 
                            s = 2.*t*u/c 

                        for k in range(i+1): 
                            t = a[k][i] 
                            u = a[k][j] 
                            a[k][i] = c*t - s*u 
                            a[k][j] = s*t + c*u
                            
                        ai = a[i]
                        aj = a[j]
                        for k in range(i+1,j): 
                            t = ai[k] 
                            u = a[k][j] 
                            ai[k] = c*t - s*u 
                            a[k][j] = s*t + c*u

                        a[j][j] = s*aij + c*ajj
                        a[i][i] = c*a[i][i] - s*(c*aij - s*ajj)

                        for k in range(j,n): 
                            t = ai[k] 
                            u = aj[k] 
                            ai[k] = c*t - s*u 
                            aj[k] = s*t + c*u
 
                        vi = v[i]
                        vj = v[j]
                        for k in range(n): 
                            t = vi[k] 
                            u = vj[k] 
                            vi[k] = c*t - s*u 
                            vj[k] = s*t + c*u 

                        a[j][i] = a[i][j] = 0.0 
                        maxd = max(maxd,abs(a[i][i]),abs(a[j][j]))
        #print iter, tol, maxdaij, nrot, '!'
        nrotsum = nrotsum + nrot
        if nrot == 0 and tol <= tolmin: 
            break
        tol = min(tol,maxdaij*1e-1,maxdaij**2)
        tol = max(tol,tolmin)

    #print "nrotsum", nrotsum
    if nrot != 0: 
        raise "Jacobi iteration did not converge in 50 passes"

    # Sort eigenvectors and values into increasing order 
    e = zerovector(n) 
    for i in range(n): 
        e[i] = a[i][i] 
        for j in range(i): 
            if e[j] > e[i]: 
                (e[i],e[j]) = (e[j],e[i]) 
                (v[i],v[j]) = (v[j],v[i]) 

    return (v,e) 

def hessian_update_bfgs(hp, dx, g, gp): 
    ''' 
    Apply the BFGS update to the approximate Hessian h[][]. 

    hp[][] = Hessian matrix from previous iteration 
    dx[]  = Step from previous iteration 
    .       (dx[] = x[] - xp[] where xp[] is the previous point) 
    g[]   = gradient at current point 
    gp[]  = gradient at previous point 

    Returns the updated hessian 
    ''' 

    n = len(hp) 
    hdx  = mxv(hp,dx) 
    dg = zerovector(n) 
    for i in range(n): 
        dg[i] = g[i] - gp[i] 
    
    dxhdx = dot(dx,hdx) 
    dxdx  = dot(dx,dx) 
    dxdg  = dot(dx,dg) 
    dgdg  = dot(dg,dg) 
    h = copymatrix(hp) 

    if (dxdx > 0.0) and (dgdg > 0.0) and (abs(dxdg/sqrt(dxdx*dgdg)) > 1.e-8): 
        for i in range(n): 
            for j in range(n): 
                h[i][j] = h[i][j] + dg[i]*dg[j]/dxdg - hdx[i]*hdx[j]/dxhdx 
    else: 
        print ' BFGS not updating dxdg (%e), dgdg (%e), dxhdx (%f), dxdx(%e)' % (dxdg, dgdg, dxhdx, dxdx) 

    return h 

def quasinr(func, guess, tol, eps, printvar=None, prnt=1): 
    ''' 
    Unconstrained minimization of a function of n variables 
    without analytic derivatives using quasi-Newtwon with BFGS update 
    and numerical gradients. 
    
    func(x) is a function that takes an array of n values and 
    returns the function value 
    
    guess[] is a vector of n values for the initial guess 
    
    tol is the convergence criterion for the maximum value 
    of the gradient 
    
    eps is the expected precision in the function value 
    
    printvar(x) is an optional user function to print the values of 
    parameters each macro iteration 
    ''' 
    
    n = len(guess) 
    x = copyvector(guess) 
    s = zerovector(n) 
    g = zerovector(n) 
    gp = zerovector(n) 
    step = zerovector(n) 
    hessian = zeromatrix(n,n) 

    alpha = 0.0 

    for iter in range(50*n): 
        (value,g,h) = numderiv(func, x, step, eps) 
        gnrm = vector_norm(g)
        
        if prnt:
            print ' ' 
            print ' iter    gnrm         value ' 
            print ' ---- --------- ----------------' 
            print "%4i %9.2e %16.8e" % (iter,gnrm,value) 

        if (prnt and printvar): 
            printvar(x) 

        if gnrm < tol: 
            if prnt: print ' Converged!' 
            break 

        if iter == 0: 
            for i in range(n): 
                hessian[i][i] = max(abs(h[i]),1e-4) 
        else: 
            hessian = hessian_update_bfgs(hessian, s, g, gp) 

        (v,e) = jacobi(hessian) 
        emax = max(map(abs,e)) 
        emin = emax*eps*100.0   # Control noise in small eigenvalues
        if prnt:
            print '\n Eigenvalues of the Hessian:' 
            printvector(e) 

        # Transform to spectral form, take step, transform back 
        gs = mxv(v,g) 
        for i in range(n): 
            if e[i] < emin:
                if prnt:
                    print ' Mode %d: small/neg eigenvalue (%f).' % (i, e[i]) 
                s[i] = -gs[i]/emin 
            else: 
                s[i] = -gs[i]/e[i] 
        s = mxv(transpose(v),s) 

        # Apply overall step restriction ... better LS obviates this 
        scale = 1.0 
        for i in range(n): 
            trust = abs(x[i])
            if trust == 0.0: trust = 0.1 # ????
            if abs(s[i]) > trust:
                if prnt:
                    print ' restricting ', i, trust, abs(x[i]), \
                          abs(x[i]/sqrt(abs(hessian[i][i]))), s[i] 
                scale = min(scale,trust/abs(s[i])) 

        if scale != 1.0: 
            for i in range(n): 
                s[i] = s[i]*scale 
        
        (alpha,value) = linesearch(func, x, s, dot(s,g), eps, prnt) 

        if alpha == 0.0:
            if prnt:
                print ' Insufficient precision to proceed further' 
            break 

        for i in range(n): 
            s[i] = s[i]*alpha 
            x[i] = x[i] + s[i] 
            gp[i]= g[i] 
    return (gnrm,x) 


def cgminold(func, dfunc, guess, tol): 
    ''' 
    Simple conjugate gradient assuming analtyic derivatives. 
    ''' 
    n = len(guess) 
    x = copyvector(guess) 
    s = zerovector(n) 
    g = zerovector(n) 
    gp= zerovector(n) 
    value = func(x) 

    for iter in range(10*n): 
        g = dfunc(x) 
        gmax = max(map(abs,g)) 
        
        print ' ' 
        print ' iter    gmax         value ' 
        print ' ---- --------- ----------------' 
        print "%4i %9.2e %16.8f" % (iter,gmax,value) 
        
        if gmax < tol: 
            print ' Converged!' 
            break 
        
        if (iter == 0) or ((iter%20) == 0): 
            beta = 0.0 
        else: 
            beta = (dot(g,g) - dot(g,gp))/(dot(s,g)-dot(s,gp)) 

        for i in range(n): 
            s[i] = -g[i] + beta*s[i] 

        (alpha,value) = linesearch(func, x, s, dot(s,g), 1e-12) 

        for i in range(n): 
            s[i] = s[i]*alpha 
            x[i] = x[i] + s[i] 
            gp[i]= g[i] 

    return (value,x) 



def cgmin(func, dfunc, guess, tol, precond=None, reset=None): 
    ''' 
    Conjugate gradient with optional preconditioning and 
    use of analytic gradients. 
    ''' 
    n = len(guess) 
    x = copyvector(guess) 
    s = zerovector(n) 
    g = zerovector(n) 
    gp= zerovector(n) 
    value = func(x) 

    if not reset: 
        reset = n 
    reset = min(reset,n) 

    for iter in range(10*n): 
        g = dfunc(x) 
        gmax = max(map(abs,g)) 
        
        print ' ' 
        print ' iter    gmax         value ' 
        print ' ---- --------- ----------------' 
        print "%4i %9.2e %16.8f" % (iter,gmax,value) 

        if gmax < tol: 
            print ' Converged!' 
            break 
        
        if precond: 
            precondg = precond(x,g) 
        else: 
            precondg = g 
        
        if (iter % reset) == 0: 
            beta = 0.0 
        else: 
            beta = (dot(precondg,g) - dot(precondg,gp))/(dot(s,g)-dot(s,gp)) 

        for i in range(n): 
            s[i] = -precondg[i] + beta*s[i] 

        (alpha,value) = linesearch(func, x, s, dot(s,g), 
                                   max(1e-16,abs(value)*1e-12)) 

        for i in range(n): 
            s[i] = s[i]*alpha 
            x[i] = x[i] + s[i] 
            gp[i]= g[i] 
        
    return (value,x) 

def cgmin2(func, guess, tol, eps, printvar=None,reset=None): 
    ''' 
    Unconstrained minimization of a function of n variables 
    without analytic derivatives using conjugate gradient with 
    diagonal preconditioning. 
    
    func(x) is a function that takes a vector of n values and 
    returns the function value 
    
    guess[] is a vector of n values for the initial guess 
    
    tol is the convergence criterion for the maximum value 
    of the gradient 
    
    eps is the expected precision in the function value 
    
    printvar(x) is an optional user function to print the values of 
    parameters each iteration 
    
    reset is the number of iterations between forced resets of the 
    conjugacy.  In principle this could be n but noise in the 
    numerical gradients makes a smaller number a better choice. 
    ''' 
    
    n = len(guess) 
    x = copyvector(guess) 
    s = zerovector(n) 
    g = zerovector(n) 
    gp = zerovector(n) 
    step = zerovector(n) 
    precondg = zerovector(n) 

    alpha = 0.0 

    if not reset: 
        reset = n 
    reset = min(reset,n) 

    for iter in range(50*n): 
        (value,g,hh) = numderiv(func, x, step, eps) 
        gmax = max(map(abs,g)) 
        
        print ' ' 
        print ' iter    gmax         value ' 
        print ' ---- --------- ----------------' 
        print "%4i %9.2e %16.8f" % (iter,gmax,value) 

        if (printvar): 
            printvar(x) 

        if gmax < tol: 
            print ' Converged!' 
            break 

        if (iter % reset) == 0: 
            # On the first iteration or if not applying conjugacy 
            # we can recompute the diagonal preconditioner 
            h = copyvector(hh) 
            for i in range(n): 
                h[i] = max(abs(h[i]),1e-6) 
        
        # Preconditioning with the diagonal of the Hessian 
        for i in range(n): 
            precondg[i] = g[i] / h[i] 
            
        # Should be able to reset every n steps but noisy gradients 
        # means that we don't have enough info. 
        if (iter % reset) == 0: 
            if iter != 0: 
                print" Resetting conjugacy" 
            beta = 0.0 
        else: 
            beta = (dot(precondg,g) - dot(precondg,gp))/(dot(s,g)-dot(s,gp)) 

        for i in range(n): 
            s[i] = -precondg[i] + beta*s[i] 

        (alpha,value) = linesearch(func, x, s, dot(s,g), eps) 

        if alpha == 0.0: 
            # LS failed, probably due to lack of precision. 
            if beta != 0.0: 
                print "LS failed - trying preconditioned steepest descent direction" 
                for i in range(n): 
                    s[i] = -g[i] 
                (alpha,value) = linesearch(func, x, s, dot(s,g), eps) 
            if alpha == 0.0: 
                print " Insufficient precision to proceed further" 
                break 

        for i in range(n): 
            s[i] = s[i]*alpha 
            x[i] = x[i] + s[i] 
            gp[i]= g[i] 
    return (value,x) 
    
def choleski(A):
    '''
    Returns the Choleski factorization of a square positive definite
    symmetric matrix.  Lij i<=j and A=L.Lt
    Raises an exception if A is not symmetric or is singular.
    '''
    n = len(A)
    m = len(A[0])
    if n != m:
        raise ValueError,"choleski factorization requires square matrix"
    for i in range(n):
        for j in range(n):
            aij = A[i][j]
            aji = A[j][i]
            if abs(aij - aji) > (1e-15*max(abs(aij),abs(aji))):
                raise ValueError,"choleski factorization requires symmetric matrix"
    L = zeromatrix(n,n)
    for k in range(n):
        sum = 0.0
        for m in range(k):
            sum = sum + L[k][m]*L[k][m]
        if A[k][k] <= sum:
            raise ValueError, "choleski factorization requires positive definite matrix"
        L[k][k] = sqrt(A[k][k]-sum)
        for i in range(k+1,n):
            sum = 0.0
            for m in range(k):
                sum = sum + L[i][m]*L[k][m]
            L[k][i] = L[i][k] = (A[i][k] - sum)/L[k][k]
    return L

def forward_elimination(L,b):
    '''
    In solving LUx=b, the first step is the forward elimination, y=L^-1.b
    '''
    n = len(L)
    y = zerovector(n)
    y[n-1] = b[n-1]/L[n-1][n-1]
    for i in range(n):
        sum = 0.0
        for j in range(i):
            sum = sum + L[i][j]*y[j]
        y[i] = (b[i] - sum)/L[i][i]
    return y

def backward_substitution(U,y):
    '''
    In solving LUx=b, the second step is the backward substitution, x=U^-1y
    '''
    n = len(U)
    x = zerovector(n)
    x[n-1] = y[n-1]/U[n-1][n-1]
    for i in range(n-2,-1,-1):
        sum = 0.0
        for j in range(i+1,n):
            sum = sum + U[i][j]*x[j]
        x[i] = (y[i] - sum)/U[i][i]
    return x

def choleski_solve(A,b):
    '''
    Solve Ax=b using Choleski factorization of A.
    '''
    L = choleski(A)
    y = forward_elimination(L,b)
    return backward_substitution(L,y)

def davidson(A,thresh=1e-6,maxsub=10,guess=None):
    '''
    Return the lowest eigenvalue and corresponding vector of A
    '''
    n = len(A)
    
    if guess:
        x = [copyvector(guess)]
        scale = 1.0/sqrt(vector_norm(x[0]))
        for i in range(n):
            x[0][i] = x[0][i]*scale
    else:
        x = [zerovector(n)]
        imin = 0
        for i in range(1,n):
            if A[i] < A[imin]:
                imin = i
        x[0][imin] = 1.0

    maxsub = min(maxsub,n)
    Ax = []
    nsub = 1
    for iter in range(10*n):
        Ax.append(mxv(A, x[nsub-1]))
        xAx = mxm(x,transpose(Ax))
        for j in range(nsub-1):
            for k in range(j+1,nsub):
                xAx[k][j] = xAx[j][k]
        #print 'Reduced matrix'
        #printmatrix(xAx)
        (v,e) = jacobi(xAx)
        #print 'Reduced evecs'
        #printmatrix(v)
        #print 'Reduced evals'
        #printvector(e)

        z = zerovector(n)
        bestx = zerovector(n)
        err = 0.0
        for i in range(n):
            Axi = 0.0
            xi  = 0.0
            for j in range(nsub):
                xi = xi + v[0][j]* x[j][i]
                Axi=Axi + v[0][j]*Ax[j][i]
            denom = (A[i][i] - e[0])
            if (denom < 0.01):
                denom = 0.01
            step = -(Axi - e[0]*xi)/denom
            #print i, step, Axi, e[0], xi, denom
            err = err + step*step
            z[i] = xi + step
            bestx[i] = xi
        scale = 1.0/sqrt(vector_norm(bestx))
        for i in range(n):
            bestx[i] = bestx[i]*scale
        print "%5d %5d %20.10f %9.1e" % (iter, nsub, e[0], sqrt(err))
        if sqrt(err) < thresh:
            return (bestx,e[0])
        
        for loop in range(2):        
            for j in range(nsub):
                zj = dot(z,x[j])
                for i in range(n):
                    z[i] = z[i] - zj*x[j][i]
        scale = 1.0/vector_norm(z)
        for i in range(n):
            z[i] = z[i]*scale
        if nsub < maxsub:
            x.append(z)
            nsub = nsub + 1
        else:
            nsub = 1
            x = [bestx]
            Ax = []
    raise "davidson did not converge"

if __name__ == '__main__':

    def precond(g): 
        # Used to test optional preconditioner for cgmin(). 
        precondg = copyvector(g) 
        for i in range(len(g)): 
            precondg[i] = precondg[i]/(i+2.0) 
        return precondg 

    def df(x): 
        d = zerovector(len(x)) 
        for i in range(len(x)): 
            d[i] = x[i]*(i+1) 
            for j in range(len(x)): 
                d[i] = d[i] + x[j] 
        return d 

    def f(x): 
        sum = 0.0 
        for i in range(len(x)): 
            for j in range(len(x)): 
                sum = sum + 0.5*x[i]*x[j] 
        for i in range(len(x)): 
            sum = sum + 0.5*x[i]*x[i]*(i+1) 
        return sum

    import random

    print '\n\n   TESTING DAVIDSON DIAGONALIZATION'
    n = 8
    a = zeromatrix(n,n)
    for i in range(n):
        for j in range(i,n):
            a[j][i] = a[i][j] = random.random()-0.5
        a[i][i] = a[i][i]*n

    (v,e) = jacobi(a)
    print ' Eval from jacobi', e[0]
    davidson(a)
    sys.exit()
    
    print '\n\n   TESTING SPARSE MATRIX PRODUCT'
    n = 13
    m = 17
    k = 11
    a = zeromatrix(k,n)
    b = zeromatrix(k,m)
    for i in range(k):
        for j in range(n):
            value = random.random()
            if (value < 0.25): a[i][j] = value
        for j in range(m):
            value = random.random()
            if (value < 0.25): b[i][j] = value
    cc= mxm(transpose(a),b)
    c = sparse_mtxm(a,b)
    err = 0.0
    for i in range(n):
        for j in range(m):
            err = err + abs(c[i][j]-cc[i][j])
    print "sparse flops:", sparse_mtxm_nflop, "dense flops:", n*m*k
    if (err < 1e-14):
        print "  OK"
    else:
        print " FAILED"
    import sys
    sys.exit()

    print '\n\n   TESTING JACOBI EIGENSOLVER\n\n' 
    n = 5
    a = zeromatrix(n,n)
    for i in range(n): 
        for j in range(i,n): 
            a[j][i] = a[i][j] = (i*j+1.)/(i+j+1.) 

    (v,e)= jacobi(a) 

    print ' eigenvalues' 
    printvector(e) 
    #print ' v ' 
    #printmatrix(v) 
    ev = mxm(v,a) 
    for i in range(n):
        norm = dot(v[i],v[i])
        if abs(norm-1.0) > 1e-14:
            print ' Error in eigenvector norm', i, norm
        etest = dot(v[i],ev[i])
        if abs(etest-e[i]) > 1e-14*max(abs(e[0]),abs(e[-1])):
            print ' Error in eigenvalue ', i, e[i], etest
        err = 0.0
        for j in range(n): 
            err = max(err,abs(ev[i][j] - e[i]*v[i][j])) 
        err = err/(n*max(1.0,abs(e[i]))) 
        if err > 1e-12: 
            print ' Error in eigenvector ', i, err
            
    sys.exit()

    print '\n\n   TESTING QUASI-NR SOLVER \n\n' 
    quasinr(f, [1.,0.5,0.3,-0.4], 1e-4, 1e-10) 

    print '\n\n   TESTING GC WITH NUM. GRAD. AND DIAG. PRECOND.\n\n' 
    cgmin2(f, [1.,0.5,0.3,-0.4], 1e-4, 1e-10, reset=20) 

    print '\n\n   TESTING GC WITH ANAL. GRAD. AND WITHOUT OPTIONAL PRECOND.\n\n' 
    cgmin(f, df, [1.,0.5,0.3,-0.4], 1e-4) 

    print '\n\n   TESTING GC WITH ANAL. GRAD. AND WITH OPTIONAL PRECOND.\n\n' 
    cgmin(f, df, [1.,0.5,0.3,-0.4], 1e-4, precond=precond) 

    print '\n\n   TESTING GC WITH ANAL. GRAD. AND NO PRECOND.\n\n' 
    cgminold(f, df, [1.,0.5,0.3,-0.4], 1e-4) 

    print '\n\n    TESTING THE CHOLESKI FACTORIZATION\n\n'
    n = 25
    A = zeromatrix(n,n)
    b = zerovector(n)
    for i in range(n):
        A[i][i] = 10.*(i+1)     # Construct A to be positive definite
        b[i] = 1.0/(i+1.0)
        for j in range(i+1,n):
            A[j][i] = A[i][j] = (i+j)/(2.0+n)
    x = choleski_solve(A,b)
    Ax = mxv(A,x)
    err = 0.0
    for i in range(n):
        err = err + abs(Ax[i] - b[i])
    print ' Cholesky linear equation solution error is ', err
