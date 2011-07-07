from math import *
from longfloat import longfloat

class Quadrature:
    def __init__(self, order, region=[-1.,1.], rule='GaussLegendre',
                 uselongfloat=0):
        self.uselongfloat = uselongfloat

        if uselongfloat:
            region = [longfloat(region[0]),longfloat(region[1])]
            
        self.order = order
        self.region = region
        self.midpoint = (region[1]+region[0])/2.
        self.scale = (region[1]-region[0])/2.
        self.rule = rule

        if rule == 'GaussLegendre':
            x,w = self.__grule()
            for i in range(self.order):
                x[i] = x[i]*self.scale + self.midpoint
                w[i] = w[i]*self.scale
        elif rule == 'GaussHermite':
            raise "sorry this is not tested"
            x,w = self.__hrule()
        elif rule == 'GaussLaguerre':
            raise "sorry this is not tested"
            x,w = self.__lrule()
        else:
            raise "unknown rule",rule

        self.x = x   
        self.w = w

    def __pn(self, n,x):
        '''
        Evaluate the Legendre polynomials up to the given order(n)
        at x in [-1,1]
        '''
        p = range(n+1)
        p[0] = 1
        if n == 0:
            return p
        p[1] = x
        for i in range(1,n):
            p[i+1] = i*(x*p[i] - p[i-1])/(i+1) + x*p[i]
        return p

    def __hn(self, n, x):
        '''
        Evaluate the renormalized Hermite polynomials Hn(x).

        Regular Hn
        p[0] = 1
        p[1] = 2*x
        p[i+1] = 2*x*p[i] - 2*i*p[i-1]

        Renormalized
        p[0] = 1/pi**(1/4)
        p[1] = sqrt(2)*x*p[0]
        p[i+1] = x*sqrt(2/(i+1))*p[i] - sqrt(j/(j+1))*p[j-1]
        
        '''
        p = range(n+1)
        p[0] = 1.0
        if self.uselongfloat:
            x = longfloat(x)
            p[0] = longfloat(1)
        if n == 0:
            return p
        p[1] = 2*x
        for i in range(1,n):
            p[i+1] = 2*x*p[i] - 2*i*p[i-1]
        return p

    def __ln(self, n, x):
        '''
        Laguerre polyn
        '''
        p = range(n+1)
        p[0] = 1.0
        if self.uselongfloat:
            x = longfloat(x)
            p[0] = longfloat(1)
        if n == 0:
            return p
        p[1] = 1-x
        for i in range(1,n):
            p[i+1] = ((2*i+1-x)*p[i]-i*p[i-1])/(i+1)
        return p
        

    def __grule(self):
        ''' Return the Gauss Legendre quadrature weights. '''

        n = self.order

        nbits = 52
        if self.uselongfloat:
            nbits = longfloat.nbits
        acc = 0.5**(nbits/3 + 10)
        if self.uselongfloat:
            acc = longfloat(acc)

        # References made to the equation numbers in Davis & Rabinowitz 2nd ed.
        roots   = range(n)
        weights = range(n)
        for k in range(n):
            # Initial guess for the root using 2.7.3.3b
            #x = (1.0 - 1.0/(8*n*n) + 1.0/(8*n*n*n))*cos(((4*k+3.0)/(4*n+2))*pi)
            x = (1.0 - 1.0/(8.0*n*n) + 1.0/(8.0*n*n*n))*cos(((4*k+3.0)/(4*n+2))*pi)
            if self.uselongfloat:
                x = longfloat(x)
            # Refine the initial guess using a 3rd-order Newton correction 2.7.3.9
            # Running this six times will ensure enough precision even if we're
            # requesting over 1500 significant decimal figures.
            for iter in range(7):
                p = self.__pn(n,x)
                p0 = p[n]                         # Value
                p1 = n*(p[n-1] - x*p[n])/(1-x*x)  # First derivative ... 2.7.3.5
                p2 = (2*x*p1 - n*(n+1)*p0)/(1-x*x)# Second derivative
                delta = - (p0/p1)*(1 + (p0*p2)/(2*p1*p1))
                x = x + delta
                if abs(delta) < acc:
                    break
            if iter >= 6:
                raise "netwon iteration for root failed"

            # Compute the weight using 2.7.3.8
            p = self.__pn(n,x)
            w = 2*(1-x*x)/(n*n*p[n-1]*p[n-1])

            roots[k] = x
            weights[k] = w

        return (roots,weights)

    def __factorial(self,n):
        f = 1.0
        if self.uselongfloat:
            f = longfloat(1.0)
        for i in range(1,n+1):
            f = f*i
        return f

    def __hrule(self):
        '''

        Gauss-Hermite rule on [0,inf]

        '''
        n = self.order

        nbits = 52
        rootpi = sqrt(pi)
        two = 2.0
        zero = 0.0
        if self.uselongfloat:
            two = longfloat(2)
            zero = longfloat(0)
            rootpi = longfloat(1).pi().sqrt()
            nbits = longfloat.nbits
        acc = 0.5**(nbits/3 + 10)
        if self.uselongfloat:
            acc = longfloat(acc)

        # Find the points which are the zeros of Hn
        # dH[n]/dx = 2*n*H[n-1]

        npt = (n-1)/2+1
        roots = range(npt)
        weights = range(npt)
        nfound = 0
        if n%2:
            p = self.__hn(n+1,0.0)
            roots[0] = zero
            weights[0] = two**(n+1)*self.__factorial(n)*rootpi/p[n+1]**2
            # !!!! Since just doing [0,inf] need to halve this
            weights[0] = weights[0]/two
            nfound = nfound + 1

        guess = 0.01
        if self.uselongfloat:
            guess = longfloat(guess)
        while nfound != npt:

            # Stupid linear search forward for sign change
            x = guess
            p = self.__hn(n,x)
            hh = p[n]
            while 1:
                x = x + 0.05
                p = self.__hn(n,x)
                if hh*p[n] < 0.0:
                    break
                
            # Cubic newton
            for iter in range(7):
                p = self.__hn(n,x)
                p0 = p[n]               # Value
                p1 = 2*n*p[n-1]         # First derivative
                p2 = 4*n*(n-1)*p[n-2]   # Second derivative
                delta = - (p0/p1)*(1 + (p0*p2)/(2*p1*p1))
                #print "newton", x, delta
                x = x + delta
                if abs(delta) < acc:
                    break
            if iter >= 6:
                raise "Hermite netwon iteration for root failed"

            new = 1
            for i in range(nfound):
                if abs(x-roots[i]) < 1e-3:
                    new = 0
                    break
            if not new:
                guess = guess + 0.1
                if guess > 100.0:
                    raise "Hermite root finder is lost"
                continue

            #print "found root", x
            roots[nfound] = x
            p = self.__hn(n+1,x)
            weights[nfound] = two**(n+1)*self.__factorial(n)*rootpi/p[n+1]**2
            nfound = nfound + 1

            guess = x + 0.05

        return roots, weights

    def __lrule(self):
        '''

        Gauss-Laguerre rule on [0,inf]

        '''
        n = self.order

        nbits = 52
        zero = 0.0
        if self.uselongfloat:
            zero = longfloat(0)
            nbits = longfloat.nbits
        acc = 0.5**(nbits/3 + 10)
        if self.uselongfloat:
            acc = longfloat(acc)

        # Find the points which are the zeros of Ln
        # dL[n]/dx   =  n*(L[n] - L[n-1])/x

        roots = range(n)
        weights = range(n)
        nfound = 0
        guess = 0.0001
        step = 0.0001
        if self.uselongfloat:
            guess = longfloat(guess)
        while nfound != n:
            # Stupid linear search forward for sign change
            x = guess
            p = self.__ln(n,x)
            hh = p[n]
            if nfound > 1:
                step = 0.1*(roots[nfound-1]-roots[nfound-2])
            print "step", step
            niter = 1000
            while niter:
                niter = niter - 1
                xp = x
                x = x + step
                p = self.__ln(n,x)
                if hh*p[n] < 0.0:
                    break
            x = (x + xp)/2.0
            if niter == 0:
                print nfound, x
                raise "Gauss Laguerre lost a root?"

            print "linear search found", x, hh, p[n]
                
            # Cubic newton
            for iter in range(7):
                p = self.__ln(n,x)
                p0 = p[n]               # Value
                p1 = n*(p[n] - p[n-1])/x # First derivative
                pp1 = (n-1)*(p[n-1] - p[n-2])/x
                p2 = n*(p1 - pp1)/x - p1/x    # Second derivative
                delta = - (p0/p1)*(1 + (p0*p2)/(2*p1*p1))
                x = x + delta
                if abs(delta) < acc:
                    break
            if iter >= 6:
                raise "Laguerre netwon iteration for root failed"
            print "found",x

            new = 1
            for i in range(nfound):
                if abs(x-roots[i]) < 1e-3:
                    new = 0
                    break
            if not new:
                guess = guess + step
                if guess > 150.0:
                    raise "Laguerre root finder is lost"
                continue

            p = self.__ln(n+1,x)
            roots[nfound] = x
            weights[nfound] = x/((n+1)*p[n+1])**2
                
            nfound = nfound + 1

            guess = x + 0.05

        return roots, weights

    def integrate(self, function=None, values=None):
        n = len(self.x)
        if function:
            values = range(n)
            for i in range(n):
                values[i] = function(self.x[i])
        sum = 0.0
        for i in range(n):
            sum = sum + self.w[i]*values[i]
        return sum

    def points(self):
        return self.x

    def weights(self):
        return self.w

def QuadratureTest(uselongfloat=0):
    ''' Test the Gauss-Legendre quadrature on [0.5,1] '''
    maxmaxerr = 0.0
    for order in range(1,20):
        q = Quadrature(order,[0.5,1.0],uselongfloat=uselongfloat)
        maxerr = 0.0
        for power in range(2*order):
            fstr = 'def __qtestf(x): return x**%d' % power
            exec(fstr)
            value = q.integrate(__qtestf)
            del(__qtestf)
            if uselongfloat:
                test = (longfloat(1) - longfloat((-power-1,1)))/(power+1.0)
            else:
                test = (1.0 - 0.5**(power+1))/(power+1.0)
            maxerr = max(maxerr,abs(float(value-test)))
        print "order %d:  integrates powers up to %d with maxerr=%e" % \
              (order, 2*order-1, maxerr)
        maxmaxerr = max(maxmaxerr,maxerr)
    print (' Maximum error from all orders %9.1e' % float(maxmaxerr))

def QuadratureTest2(uselongfloat=0):
    ''' Test the Gauss-Hermite quadrature on [-inf..inf] '''
    maxmaxerr = 0.0
    rootpi = sqrt(pi)
    two = 2.0
    if uselongfloat:
        two = longfloat(2)
        rootpi = two.pi().sqrt()

    for order in range(2,24):
        q = Quadrature(order,rule="GaussHermite",
                       uselongfloat=uselongfloat)
        maxerr = 0.0
        for power in range(0,2*order,2):
            fstr = 'def __qtestf(x): return x**%d' % power
            exec(fstr)
            value = q.integrate(__qtestf)
            del(__qtestf)
            test = 0.5*rootpi*0.5**(power/2)
            for i in range(power-1,0,-2):
                test = test * i
            maxerr = max(maxerr,abs(float(value-test)/test))
        print "order %d:  integrates powers up to %d with maxerr=%e" % \
              (order, 2*order-2, maxerr)
        maxmaxerr = max(maxmaxerr,maxerr)
    print (' Maximum error from all orders %9.1e' % float(maxmaxerr))


if __name__ == '__main__':

    q = Quadrature(7,(0.0,1.0),rule='GaussLegendre')
    print q.points()
    print q.weights()
    stop

    print ' \n\n Testing quadrature with default floats \n'
    QuadratureTest()

    print ' \n\n Testing quadrature with long floats \n'
    QuadratureTest(uselongfloat=1)
    
    #QuadratureTest2()

    #QuadratureTest2(1)
