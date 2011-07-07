from math import *
from types import *
import re

# Outstanding issues 
#   Convergence criteria for series expansion of functions
#   How to round correctly?
#   Also, rounding while printing
#   Faster truncation
#   Use of interp routines for function evaluation

class interptable:
    def __init__(self, func, lo, hi, n, precomputed=None):
        '''
        Tablulate the function over n intervals or n+1 points
        in [lo,hi] inclusive of the end points.

        This is now hardwired to longfloat, but the changes
        to move it back to generality are obvious.
        '''
        if precomputed:
            (lo, hi, n, h, h2, scale, func, x, f, nbits) = precomputed
        else:
            lo = longfloat(lo)
            hi = longfloat(hi)
            h = (hi - lo)/(n+1)
            h2 = h*0.5
            scale = 1.0/h
            x = [0]*(n+1)
            f = [0]*(n+1)
            for i in range(n+1):
                x[i] = lo + i*h
                f[i] = func(x[i])
            nbits = longfloat.nbits
        self.lo = lo
        self.hi = hi
        self.n = n
        self.h = h
        self.h2 = h2
        self.scale = scale
        self.func = func
        self.x = x
        self.f = f
        self.nbits = nbits

    def map(self,x):
        '''
        If x is within the range of the interpolation table return
        k - the index of the nearest node x[k]
        d - the difference x-x[k]

        If x is out of range raise a value exception.
        '''
        if self.lo <= x <= self.hi:
            k = int((x-self.lo)*self.scale)
            d = x - self.x[k]
            if d > self.h2:
                k = k + 1
                d = x - self.x[k]
            return (k,d)
        else:
            raise ValueError, ('when interpolating %s' % self.func)

    def __repr__(self):
        '''
        This so that tables may be precomputed for fast startup.
        '''
        return 'interptable(None,None,None,None,precomputed=' + \
              repr((self.lo,self.hi,self.n,self.h,self.h2,self.scale,
                    str(self.func),self.x,self.f,self.nbits)) + ')'

class interpolators:
    def __init__(self):
        # Exp tabulated on [-1/2,1/2]
        half = longfloat((-1,1))
        # 20 .12
        # 40 .12
        #100 .10
        #400 .088
        ###self.texp = interptable(lambda (x):x.exp(), -half, half, 40)
        #print self.texp

        # Sin and Cos on [0,pi/2] (1.571 is slightly bigger than pi/2)
        zero = longfloat(0)
        pi2  = longfloat(1.571)
        ###self.tsin = interptable(lambda (x):x.sin(), zero, pi2, 30)
        ###self.tcos = interptable(lambda (x):x.cos(), zero, pi2, 30)

        # log(x) with x in [1/sqrt(2), sqrt(2)]
        # 1.415 > sqrt(2)    0.707 < 1/sqrt(2)
        n = 100
        self.tlog = interptable(lambda (x):x.log(),
                                longfloat(0.707), longfloat(1.415), n)
        self.tlog.rx = [0]*(n+1)
        for k in range(n+1):
            self.tlog.rx[k] = self.tlog.x[k].reciprocal()

    def exp(self,x):
        '''
        For arguments in the range [0,1/2] exp(x) is O(1) so it
        suffices to ensure that the remainder is < eps
        '''
        (k, d) = self.texp.map(x)
        #print k , d
        fk = self.texp.f[k]
        result = fk
        term = fk
        eps = x.eps()
        for i in range(1,50):
            term = term * d / i
            result = result + term
            if abs(term) < eps:
                #print ' exp terminating at ', i
                return result
        raise ArithmeticError, "interp exp did not converge"

    def sin(self,x):
        '''
        For arguments in the range [0,pi/2] sin(x) is O(min(x,1)).
        So, unless x is a hard zero, use relative error.
        '''
        if not x:
            return longfloat(0)
        (k, d) = self.tsin.map(x)
        sink = self.tsin.f[k]
        cosk = self.tcos.f[k]
        result = sink
        term = d
        eps = x.eps()
        if x < 1:
            eps = eps*x
        for i in range(1,50,4):
            result = result + cosk*term
            term = term*d/(i+1)
            result = result - sink*term
            term = term*d/(i+2)
            result = result - cosk*term
            term = term*d/(i+3)
            result = result + sink*term
            term = term*d/(i+4)
            if abs(term) < eps:
                #print ' sin terminating at ', i+3
                return result
        raise ArithmeticError, "interp sin did not converge"
        
    def cos(self,x):
        '''
        Like sin, use relative error.
        '''
        if not x:
            return longfloat(1)
        (k, d) = self.tcos.map(x)
        sink = self.tsin.f[k]
        cosk = self.tcos.f[k]
        result = cosk
        term = d
        eps = x.eps()
        if x < 1:
            eps = eps*x
        for i in range(1,50,4):
            result = result - sink*term
            term = term*d/(i+1)
            result = result - cosk*term
            term = term*d/(i+2)
            result = result + sink*term
            term = term*d/(i+3)
            result = result + cosk*term
            term = term*d/(i+4)
            if abs(term) < eps:
                #print ' cos terminating at ', i+3
                return result
        raise ArithmeticError, "interp sin did not converge"

    def log(self, x):
        '''
        log(x+d) = log(x) - sum(i=1..inf) (-d/x)^i/i

        For x in the range [1/root(2),root(2)] log(x) is about
        O(x-1) ... so use relative error to this for convergence.

        log(x) = 2*(p + p^3/3 + p^5/5 + ...) where p = (x-1)/(x+1)

        log(xk+d) = log(xk)+log(1+d/xk)

        '''
        if x == 1:
            return longfloat(0)
        (k, d) = self.tlog.map(x)
        dx = d*self.tlog.rx[k]
#        print k, float(d), float(self.tlog.x[k]), float(dx)
        p = dx/(dx+longfloat(2))
#        print ' p', float(p)
        p2 = p*p
        eps = x.eps() * abs(self.tlog.f[k])
#        print ' eps ', float(eps)
        term = longfloat(1)
        result = longfloat(0)
        for i in range(1,501,2):
            result = result + term/longfloat(i)
            term = term*p2
            if abs(term) < eps:
                break
        if i >= 498:
            raise ArithmeticError, "interp log did not converge"
        return (2*p*result + self.tlog.f[k])
        
##         dx = d/x
##         result = self.tlog.f[k]
##         term = dx
##         eps = x.eps() * abs(x-1) # relative error
##         for i in range(1,50):
##             result = result + term/i
##             term = term * dx
##             if abs(term) < eps:
##                 #print ' log terminating at ', i
##                 return result
##         raise ArithmeticError, "interp log did not converge"

    def sindsin(self, x):
        return (self.sin(x), self.cos(x))

    def cosdcos(self, x):
        return (self.cos(x), -self.sin(x))

    def expdexp(self, x):
        return (self.exp(x), self.exp(x))

    def logdlog(self, x):
        return (self.log(x), 1/x)

    def invert(self, y, gdg, x):
        '''
        solve for x in y = g(x) ... i.e., invert g assuming that
        we can compute both g and its derivative.

        x inputs an initial guess for x which had better be in the
        quadratically convergent region.
        '''
        # delta = (g(x) - y)/g'(x)
        eps = x.eps()
        while 1:
            (g,dg) = gdg(x)
            print ' g dg ', g, dg
            delta = (g - y)/dg
            print ' delta ', x, delta
            x = x - delta
            if abs(delta) < eps:
                break
        return x

class longfloat:
    '''
    Represents a floating point number to arbitrary precision
    by storing the mantissa as a long integer and the exponent
    as an integer.

    value = 2**exponent * mantissa

    The class member nbits is used to control the number of bits
    of precision retained by operations.  It defaults to 208 which
    is four times that of the standard double precision (52 bits).
    It may be modified as follows
    
    longfloat.nbits = newvalue
    
    and the modified value will be adopted by all new instances.
    
    The class defines the standard methods to support
    -  conversion to/from integers, floats, longs, strings, 2-tuple
    -  addition, subtraction, multiplication, negation, division
    -  comparison and non-zero test
    -  printing and representation as a string (__repr__)

    and also defines methods to compute the following
    
    -  reciprocal, reciprocal square root and square root

    After each operation the result is truncated with rounding to
    the required finite precision.  Internally, some operations
    will be performed to a higher precision.

    The method eps() returns the smallest number such that 1+eps > 1.

    Some rather limited testing suggests that the errors are
    -  the last bit in the basic operations +*- due to rounding
    -  the last two bits in divide and square root
    -  the last two or three bits in log/exp/sin/cos
    -  the last four or five bits in asin
    -  the last five to eight bits in acos
    '''

    nbits   = 208
    debug   = 0

    __name__ = 'longfloat'

    __small_int_recip = [(0,0)]*501  # Cache 1/n  for 0<n<=500 ... entry is (1/n,nbits)
   
    def __init__(self, value, nbits=None):
        # Constructor allows override of class default for nbits
        if nbits != None:
            self.nbits = nbits
        if self.nbits <= 0:
            raise ValueError, "longfloat.nbits must be > 0"
            
        # Initialization from ints and longfloat done here, others done
        # thru coerce
        t = type(value)
        if isinstance(value,longfloat):
            self.exponent = value.exponent
            self.mantissa = value.mantissa
        elif t == IntType or t == LongType:
            self.exponent = 0
            self.mantissa = long(value)
        else:
            (junk, a) = self.__coerce__(value)
            self.exponent = a.exponent
            self.mantissa = a.mantissa

    def __power2(self, exponent):
        ' return 2.0**exponent '
        if exponent < 0:
            fac = 0.5
        else:
            fac = 2.0
        
        scale = 1.0
        for i in range(abs(exponent)):
            scale = scale * fac
        return scale
    
    def __coerce__(self, other):
        ' coerce strings, ints, longs, floats and longfloats into longfloat'
        if isinstance(other,longfloat):
            return (self, other)

        t = type(other)
        if self.debug: print ' longfloat: coercing ', t, other

        value = longfloat(0)
        if (t == IntType) or (t == LongType):
            value.exponent = 0
            value.mantissa = long(other)
        elif t == StringType:
            # The exact representation from __repr__() is a tuple of 
            # two integers (exponent, mantissa).
            # Alternatively match a floating point number
            # Match two signed integers from a 2-tuple, e.g., (1,2L)
            pat = r'\([ ]*([+-]?\d+)[lL]?[ ]*,[ ]*([+-]?\d+)[lL]?[ ]*\)'
            m = re.search(pat,other)
            if m:
                (e,m) = m.groups()
                value.exponent = int(e)
                value.mantissa = long(m)
            else:
                patf = r'\A\s*([+-]?\d+){1,1}(\.\d*)?([eE][+-]?\d+)?'
                m = re.search(patf,other)
                mantissa = ''
                if m:
                    mantissa = ''
                    e10 = 0
                    (before,after,expnt) = m.groups()
                    print ' before, after, expnt ', before, after, expnt
                    if before:
                        mantissa = mantissa + before
                    if after:
                        after = after[1:]        # Get rid of leading .
                        e10 = -len(after)
                        mantissa = mantissa + after
                    if expnt:
                        e10 = e10 + int(expnt[1:])
                    print ' mantissa ', mantissa, ' e10 ', e10
                    if e10:
                        longfloat.nbits = longfloat.nbits + 31
                        log10  = self.__log10()
                        log2   = self.__log2()
                        log210 = log10/log2
                        e2 = e10 * log210        # Base two exponent
                        r2 = e2 - int(e2)        # Fractional part
                        print ' r2 ', float(r2)
                        e2 = int(e2)             # Integer part
                        r2 = (r2 * log2).exp()   # 2^r2
                        print ' 2^r2 ', float(r2)
                        print ' r2* ', float(r2 - r2.log()/log2)
                        value.exponent = e2
                        value.mantissa = long(mantissa)
                        value = value * r2
                        longfloat.nbits = longfloat.nbits - 31
                        value.truncate()
                    else:
                        value.exponent = 0
                        value.mantissa = m
                else:
                    raise ("longfloat: fromstr: '%s' is neithr a tuple (exponent,mantissa) nor a floating point number" % other)
        elif (t == TupleType) and (len(other) == 2):
            value.exponent = int(other[0])
            value.mantissa = long(other[1])
        elif t == FloatType:
            # A double has a 52 bit mantissa and we need to preserve
            # all of those bits ... but storing too many is inefficient.
            # Must also watch out for very small numbers that will
            # overflow when reciprocated

            # This routine is paranoid and always checks that an exact
            # representation of the input was generated.

            if other == 0.0:
                exponent = 0
                mantissa = 0L
            else:
                exponent = int(log(abs(other))/log(2.0))
                exponent = exponent - 53
                extra = 0
                while exponent < -990:
                    other = other * 2.0
                    exponent = exponent + 1
                    extra = extra + 1
                mantissa = long(other*self.__power2(-exponent))
                exponent = exponent - extra
                while extra > 0:
                    other = other * 0.5
                    extra = extra - 1

                # Eliminate trailing zeros in the mantissa for
                # values that are actually exactly representable
                while (mantissa & 1L) == 0:
                    mantissa = mantissa >> 1
                    exponent = exponent + 1
                
            test = float(mantissa)*self.__power2(exponent)
            if (test - other) != 0.0:
                print ' exponent ', exponent
                print ' mantissa ', mantissa
                print ' input      %24.17e' % other
                print ' test       %24.17e' % test
                print ' err        %24.17e' % (test-other)
                raise TypeError, "failed to convert double to longfloat"
                
            value.mantissa = mantissa
            value.exponent = exponent
        else:
            print ' unknown TYPE IS ', t
            raise TypeError, 'conversion to longfloat failed'

        return (self,value)

    def __float__(self):
        # If the mantissa is too long then the standard conversion
        # of long to float does not work.  Ugh.

        e = self.exponent
        m = self.mantissa
        while abs(m) > 1e17:
            m = m >> 1
            e = e + 1
        return (float(m)*self.__power2(e))  # THIS CAN FAIL .. FIX IT!!!!

    def __int__(self):
        return int(long(self))

    def __long__(self):
        e = self.exponent
        m = self.mantissa
        if e > 0:
            m = m << e
        else:
            m = m >> abs(e)
        return m
    
    def __repr__(self):
        ' exact string representation is a tuple (exponent,mantissa) '
        return "longfloat((" + repr(self.exponent) + "," + repr(self.mantissa) + "))"

    def __str__(self):
        ' human readable but inexact string representation is decimal float'
        if not self:
            return '0'
            
        # Need to get the exponent into base 10 and scale the
        # mantissa by its fractional part.
        # 2^e * m = 10^(e/log2(10)) * m = 10^e10 * (10^r10 * m)
        # base 10 exponent is e10+r10 (integer part + remainder)

        nbits   = self.nbits
        ndigits = int((self.nbits-1)*0.30103) + 1
        longfloat.nbits = longfloat.nbits + 31 # Should be log2(|e+nbm|)?

        x = longfloat(self)
        e = x.exponent
        m = x.mantissa
        nbm = x.__nbits_in_mantissa()
        if nbm < nbits:
            e = e - (nbits - nbm)
            m = m << (nbits - nbm)
            x.exponent = e
            x.mantissa = m

        log10  = self.__log10()
        log2   = self.__log2()
        log210 = log10/log2

        # Split exponent into fractional and integer parts
        e10  = e / log210
        r10 = e10 - int(e10)
        e10 = int(e10)
        tenr10 = (r10*log10).exp()
        
        # Generate 10^r10 * m truncated to the appropriate no. of bits
        tmp = longfloat(x)
        tmp.exponent = 0
        tmp = tmp*tenr10

        # Shift the mantissa to zero the remaining exponent
        etmp = tmp.exponent
        mtmp = tmp.mantissa
        m10 = tmp.mantissa

        if etmp > 0:
            m10 = m10 << etmp
        elif etmp < 0:
            # Cannot just shift ... must round correctly.
            # Examine the last bit discarded
            etmp = abs(etmp)
            m10 = m10 >> (etmp-1)
            if (m10 & 1):
                m10 = m10 + 2
            m10 = m10 >> 1
            
        # Now we have e10 and m10
        # Count the no. of decimal figures in m10 so that we can
        # have one significant figure before the decimal point
        strm = str(abs(m10))
        n10 = len(strm) - 2
        e10 = e10 + n10

        # Finally can generate the result
        signm = '+'
        if m10 < 0: signm = '-'
        signe = '+'
        if e10 < 0: signe = '-'

        lead = strm[0]
        strm = strm[1:]
        if strm[-1] == 'L':
            strm = strm[:-1]
        while len(strm) and (strm[-1] == '0'):
            strm = strm[:-1]
        if strm == '':
            strm = '0'
        
        longfloat.nbits = longfloat.nbits - 31

        return signm + lead + '.' + strm + 'e' + signe + str(abs(e10))
        

    def __mul__(self,other):
        value = longfloat(0)
        value.mantissa = self.mantissa*other.mantissa
        value.exponent = self.exponent+other.exponent
        value.truncate()
        return value

    def __rmul__(self,other):
        return self*other

    def __common_exponent(self,other):
        '''
        Return the mantissas scaled to a common exponent
        '''
        eself = self.exponent
        mself = self.mantissa
        eother = other.exponent
        mother = other.mantissa

        if eother < eself:
            exponent = eother
            mself = mself << (eself-eother)
        else:
            exponent = eself
            mother = mother << (eother-eself)
        return (exponent, mself, mother)

    def __add__(self,other):
        (exponent, mself, mother) = self.__common_exponent(other)
        value = longfloat(0)
        value.exponent = exponent
        value.mantissa = mself + mother
        value.truncate()
        return value
            
    def __radd__(self,other):
        return self+other

    def __sub__(self,other):
        return self+(-other)

    def __rsub__(self,other):
        return other+(-self)

    def __neg__(self):
        value = longfloat(0)
        value.exponent =  self.exponent
        value.mantissa = -self.mantissa
        return value

    def __abs__(self):
        value = longfloat(0)
        value.exponent =  self.exponent
        value.mantissa = abs(self.mantissa)
        return value
        
    def __pos__(self):
        return self

    def __nonzero__(self):
        return (self.mantissa != 0L)

    def __cmp__(self, other):
        (exponent, mself, mother) = self.__common_exponent(other)
        if mself < mother:
            return -1
        elif mself == mother:
            return 0
        else:
            return 1

    def __pow__(self,power):
        if power == 0:
            return longfloat(1)
        elif power == 1:
            return longfloat(self)
        elif power == 2:
            return self*self
        elif power == 3:
            return self*self*self
        elif power == 4:
            x2 = self*self
            return x2*x2

        if power < 0 and self == 0:
            raise ValueError, "longfloat: pow: negative power of zero"

        x = longfloat(self)
        phase = 1
        if power.exponent == 0:
            # Integer power ... can handle negative arguments
            if x < 0:
                if (power % 2):
                    phase = -1
                x = -x
            # Repeated squaring for integer powers
            m = int(power.mantissa)
            if m < 0:
                x = x.reciprocal()
                m = -m
            xk = x                  # x^(2^k)
            result = 1.0
            while m:
                if (m&1):
                    result = result * xk
                xk = xk*xk
                m = m >> 1
            if phase < 0:
                result = -result
            return result
        
        if x < 0:
            raise ValueError, "longfloat: pow: negative number to a non-integer power"

        logx = x.log()
        result = (power*logx).exp()
        if phase < 0:
            result = -result

        return result

    def __rpow__(self,power):
        return power.__pow__(self)

    def rsqrt(self):
        '''
        Compute the reciprocal squareroot using Newton-Raphson iteration.
        '''
        if self <= 0.0:
            raise ValueError, 'reciprocal squareroot of non-positive value'

        nbits = self.nbits

        # Generate the initial guess using the standard float operation
        # but we need to make sure it is representable.
        nbm = self.__nbits_in_mantissa()
        e = self.exponent
        if ((e+nbm)%2):
            nbm = nbm + 1
        self.exponent = -nbm
        s = longfloat(1.0/sqrt(float(self)))
        self.exponent = e
        s.exponent = s.exponent - ((e + nbm)>>1)

        threehalf = longfloat(3)
        threehalf.exponent = threehalf.exponent - 1
        xhalf = longfloat(self)
        xhalf.exponent = xhalf.exponent - 1
        nbitsgot = 51
        while nbitsgot < nbits:
            #snew = half*s*(three - s*s*self)
            snew = s*(threehalf - s*s*xhalf)
            s = snew
            nbitsgot = nbitsgot*2-1 # Quad. conv. with rounding loss

        return s

    def sqrt(self):
        if self.mantissa == 0L:
            return longfloat(0)
        else:
            return self*self.rsqrt()

    def reciprocal(self):
        '''
        Compute the reciprocal using Newton-Raphson iteration.
        '''
        if not self:
            raise ValueError, 'reciprocal of zero value'

        if not self.exponent:
            # We have an integer.  The reciprocals of small integers
            # are precomputed and cached for efficiency
            r = self.__small_integer_reciprocal()
            if r:
                return r
        
        nbits = self.nbits

        # Generate the initial guess using the standard float operation
        # but we need to make sure it is representable.
        nbm = self.__nbits_in_mantissa()
        e = self.exponent
        self.exponent = -nbm
        r = longfloat(1.0/float(self))
        self.exponent = e
        r.exponent = r.exponent - (e + nbm)
        
        two = longfloat(2)
        nbitsgot = 51
        while nbitsgot < nbits:
            rnew = r*(two - r*self)
            r = rnew
            nbitsgot = nbitsgot*2 - 1 # Quad. conv. with rounding loss

        return r

    def __small_integer_reciprocal(self):
        ''' Maintain a cache of the reciprocals of small integers '''
        
        m = self.mantissa
        absm = abs(m)
        if absm > 500:
            return None
        absm = int(absm)
        
        (val,nbits) = longfloat.__small_int_recip[absm]
        if nbits < longfloat.nbits:
            val = 0
            
        if not val:
            #print ' computing small int ', absm
            # Not in cache. Set expnt = -1 to stop reciprocal calling
            # here in infinite loop.
            val = longfloat((-1,absm)).reciprocal()
            val.exponent = val.exponent - 1
            longfloat.__small_int_recip[absm] = (val,longfloat.nbits)
        
        # Finally, return a copy of the cached object
        val = longfloat(val).truncate()
        if m < 0:
            return -val
        else:
            return val
        
    def __log10(self):
        # Cache the value of log10 in the class but we
        # need to keep track of the precision used to compute it.
        try:
            if longfloat.__log10nbits >= self.nbits:
                return longfloat(longfloat.__log10value)
        except AttributeError: 
            pass

        longfloat.__log10value = longfloat(10).log()
        longfloat.__log10nbits = self.nbits
        return longfloat(longfloat.__log10value)

    def __log2(self):
        # Cache the value of log2 in the class but we
        # need to keep track of the precision used to compute it.
        try:
            if longfloat.__log2nbits >= self.nbits:
                return longfloat(longfloat.__log2value)
        except AttributeError: 
            pass
        
        longfloat.nbits = longfloat.nbits + 8

        test = longfloat((-self.nbits-4,1))
        x = longfloat(2)
        p = (x-1)/(x+1)
        p2 = p*p
        term = longfloat(1)
        result = longfloat(0)
        i = 1
        while 1:
            result = result + term/i
            term = term*p2
            i = i + 2
            if abs(term) < test:
                break
        log2 = longfloat.__log2value = 2*p*result
        longfloat.__log2nbits = self.nbits

        longfloat.nbits = longfloat.nbits - 8

        return longfloat(log2)

    def log(self):
        '''
        Return the natural log of self. self is not modified.

        Use log(x) = log(x * 2^-k * 2^k) = k*log(2) + log(x*2^-k)
        to scale x into (1/2) <= x <= 1.

        Then for the scaled variable use

        log(x) = 2*(p + p^3/3 + p^5/5 + ...) where p = (x-1)/(x+1)
        '''

        # CANNOT CALL PRINT OF LONGFLOAT FROM THIS ROUTINE SINCE PRINT
        # CALLS THIS !!!!!
        
        if self <= 0:
            raise ValueError, "longfloat: taking log of non-positive value"

        longfloat.nbits = longfloat.nbits + 8

        e = self.exponent
        m = self.mantissa

        # Initially scale the argument into [0.5,1], then into
        # [1/sqrt(2), sqrt(2)] which gives the smallest range for p.
        k = self.__nbits_in_mantissa()
        power = e + k

        x = longfloat((-k,m))
        x707 = longfloat(181)
        x707.exponent = -8    # 181/256 = 0.707
        if abs(x) < x707:
            x.exponent = x.exponent + 1
            power = power - 1
            
        p = (x-1)/(x+1)
       
        # Apply the standard power-series
        test = longfloat((-self.nbits-4,1))*abs(x-1)
        term = longfloat(1)
        p2 = p*p
        result = longfloat(0)
        for i in range(1,501,2):
            result = result + term/i
            term = term*p2
            if abs(term) < test:
                break
        if i >= 498:
            raise ArithmeticError, "longfloat log did not converge"

        result = 2*p*result

        longfloat.nbits = longfloat.nbits - 8

        result = power*self.__log2() + result
        result.truncate()
        return result

    def exp(self):
        '''
        Return exp(self).  self is not modified.
        '''
        # CANNOT CALL PRINT OF LONGFLOAT FROM THIS ROUTINE SINCE PRINT
        # CALLS THIS !!!!!

        if not self:
            return longfloat(1)
        
        longfloat.nbits = longfloat.nbits + 8

        e = self.exponent
        m = self.mantissa

        # Scale the argument so that it is less than 1/2 (try 1/4)
        k = self.__nbits_in_mantissa()
        k = k + 1                       # to get 1/4 rather than 1/2
        power = e + k
        if power < 0:
            power = 0

        x = longfloat(0)
        x.exponent = e - power
        x.mantissa = m
            
        # Apply the standard power-series
        test = longfloat((-self.nbits-4,1))*abs(x)
        term = longfloat(1)
        result = longfloat(0)
        for i in range(1,500):
            result = result + term
            term = term * x / i
            if abs(term) < test:
                break
        if i >= 499:
            raise ArithmeticError, "longfloat exp did not converge"
        
        while power > 0:
            result = result * result
            power = power - 1

        longfloat.nbits = longfloat.nbits - 8
        result.truncate()

        return result

    def pi(self):
        '''
        Compute pi with the quadratically convergent Borwein algorithm
        which has monotonically decreasing error pin-pi < 10^-2^(n+1)
        '''
        try:
            if longfloat.__pinbits >= self.nbits:
                p = longfloat(longfloat.__pi)
                p.truncate()
                return p
        except AttributeError: 
            pass
            
        longfloat.nbits = longfloat.nbits + 8
        
        half = longfloat(0.5)
        x = longfloat(2).sqrt()         # x0
        mypi = 2 + x                      # pi0
        y = x.sqrt()                    # y1
        rsqrtx = x.rsqrt()
        sqrtx  = x*rsqrtx
        x = half*(sqrtx + rsqrtx)    # x1

        test = longfloat((-self.nbits/2-1,1))
        for iter in range(30):
            ryn1 = (y+1).reciprocal()
            mypi_new = mypi*(x+1)*ryn1
            err = mypi - mypi_new
            mypi = mypi_new
            if err < test:
                break
            rsqrtx = x.rsqrt()
            sqrtx  = x*rsqrtx
            y = (y*sqrtx + rsqrtx)*ryn1
            x = half*(sqrtx + rsqrtx)    # x1
        if iter >= 29:
            raise ArithmeticError, 'longfloat pi did not converge'

        mypi.nbits = longfloat.nbits = longfloat.nbits - 8
        mypi.truncate()

        longfloat.__pi = mypi
        longfloat.__pinbits = mypi.nbits

        return longfloat(mypi)

    def __sin(self):
        ''' Compute sine of argument already scaled to [0,pi/4] '''

        x  = self
        x2 = x*x
        
        # Apply the standard power-series
        test = longfloat((-self.nbits-4,1))
        term = longfloat(x)
        result = longfloat(0)
        i = 1
        for i in range(1,301,2):
            result = result + term
            term = -term * x2 / (i+1) / (i+2) # Recips of small ints are fast
            if abs(term) < test:
                break
        if i >= 299:
            raise ArithmeticError, "longfloat sin did not converge"
            
        return result

    def __cos(self):
        ''' Compute cosine of argument already scaled to [0,pi/4] '''

        x  = self
        x2 = x*x
        
        # Apply the standard power-series
        test = longfloat((-self.nbits-4,1))
        term = longfloat(1)
        result = longfloat(0)
        for i in range(0,300,2):
            result = result + term
            term = -term * x2 / (i+1) / (i+2)  # Recips of small ints are fast
            if abs(term) < test:
                break
        if i >= 298:
            raise ArithmeticError, "longfloat cos did not converge"
            
        return result

    def sin(self):
        longfloat.nbits = longfloat.nbits + 8
        x = longfloat(self)
        # Shift argument into [0,2pi]
        pi    = self.pi()
        twopi = pi+pi
        pi2   = pi*0.5
        pi4   = pi*0.25
        n     = int(x/twopi)
        x     = x - n*twopi

        # Shift argument into [0,pi/2] using
        # sin(x) = sin(pi-x)=-sin(pi+x)=-sin(2pi-x)
        phase = 1
        if x > pi:
            phase = -1
            x = x - pi
        if x > pi2:
            x = pi - x

        # Shift into [0,pi/4] using
        # sin(x+pi/4) = (sin(x) + cos(x)) / sqrt(2)

        if x >= pi4:
            x = x - pi4
            result = phase * (x.__sin() + x.__cos()) * longfloat(2).rsqrt()
        else:
            result = phase * x.__sin()

        longfloat.nbits = longfloat.nbits - 8

        result.truncate()
        return result

    def cos(self):
        longfloat.nbits = longfloat.nbits + 8
        x = longfloat(self)
        # Shift argument into [0,2pi]
        pi    = self.pi()
        twopi = pi+pi
        pi2   = pi*0.5
        pi4   = pi*0.25
        n     = int(x/twopi)
        x     = x - n*twopi

        # Shift argument into [0,pi/2] using
        # cos(x) = -cos(pi-x) = -cos(pi+x)= cos(2pi-x)
        phase = 1
        if x > pi:
            x = x - pi
        if x > pi2:
            phase = -1
            x = pi - x

        # Shift into [0,pi/4] using
        # cos(x+pi/4) = (cos(x) - sin(x)) / sqrt(2)

        if x >= pi4:
            x = x - pi4
            result = phase * (x.__cos() - x.__sin()) * longfloat(2).rsqrt()
            ##cosx = x.__cos()
            ##sinx = ((1+cosx).sqrt())*((1-cosx).sqrt())
            ##result = phase * (cosx - sinx) * longfloat(2).rsqrt()
        else:
            result = phase * x.__cos()

        longfloat.nbits = longfloat.nbits - 8
        result.truncate()
        return result

    def acos(self):
        '''
        pi/2 - asin(x) -1 < x < 1
        '''
        pi  = self.pi()
        pi2 = 0.5*pi

        # Must first shift the argument into a range where the subtraction
        # pi/2 - asin(x) does not lose precision

        x = longfloat(self)
        xold = longfloat(x)
        # Restrict to 0 <= x <= 1 with
        # acos(x) = pi - acos(-x)
        extra = 0
        phase = 1
        if x < 0:
            extra = pi
            phase = -1
            x = -x

        xx = x
        # Restrict to 0 <= x <= 1/sqrt(2) with
        # acos(x) = pi/4  + acos(-x + sqrt(1-x^2))
        extra2 = 0
        rroot2 = longfloat(2).rsqrt()
        if x > rroot2:
            extra2 = 0.25*pi
            #x = (x - ((1+x).sqrt())*((1-x).sqrt()))*rroot2
            x = (x - (1-x*x).sqrt())*rroot2

        asinx = x.asin()
        #print ' asin arg ', float(x), float(asinx) - asin(float(x))
        acosx = pi2 - asinx
        #print ' acos arg ', float(x), float(acosx) - acos(float(x))
        if extra2:
            acosx = acosx - extra2
            #print ' acos xx ', float(xx), float(acosx) - acos(float(xx))
        
        acosx = extra + phase*acosx
        #print ' final ', float(xold), float(acosx) - acos(float(xold))

        return acosx

    def asin(self):
        '''
        asin(x) = sum(n) (-1)^n C(-1/2,n) x^(2n+1)/(2n+1)
        %       = x + (1/2)*x^3/3 + (1*3/2*4)*x^5/5 + (1*3*5/2*4*6)*x^7/7 ...

        asin(x) = -asin(-x)
        '''
        longfloat.nbits = longfloat.nbits + 8

        x = longfloat(self)
        one = longfloat(1)
        if (x < -one) or (x > one):
            longfloat.nbits = longfloat.nbits - 8
            raise ValueError, 'longfloat asin argument not in [-1,1]'

        input_phase = 1
        if x < 0:
            input_phase = -1
            x = -x

        pi  = self.pi()
        pi2 = pi*0.5
        if x == one:
            longfloat.nbits = longfloat.nbits - 8
            return input_phase*pi2

        # Shift the argument into [0,1/sqrt(2)] using
        # asin(x) = pi/2 - asin(sqrt(1-x^2))  0<=x<=1
        # asin(x) =-pi/2 + asin(sqrt(1-x^2)) -1<=x<=0
        # We need y = sqrt(1-x^2) accurate enough so that asin(y)
        # is accurate compared to pi/2.

        extra = 0
        phase = 1
        if x > longfloat(2).rsqrt():
            extra = pi2
            phase = -1
            #x = ((1+x).sqrt())*((1-x).sqrt())
            x = (1-x*x).sqrt()

        # Shift the argument into [0,1/2] using
        # asin(x) = pi/6 + asin(x*sqrt(3)/2 - sqrt(1-x^2)/2)

        didpi6 = 0
        half = longfloat((-1,1))
        if x > half:
            didpi6 = 1
            #x = (x*longfloat(3).sqrt() - ((1+x).sqrt())*((1-x).sqrt())) * half
            x = (x*longfloat(3).sqrt() - (1-x*x).sqrt()) * half
            #print ' shifting x to ', x, half
                 
        # asin(x) = sum(n) (-1)^n C(-1/2,n) x^(2n+1)/(2n+1)
        # %       = x + (1/2)*x^3/3 + (1*3/2*4)*x^5/5 + (1*3*5/2*4*6)*x^7/7 ...

        x2 = x*x
        term = x
        result = longfloat(0)
        test = longfloat((-self.nbits,1))

        # 0 <= x2 <= 0.25 so we get at least 2 bits per iteration
        maxi = self.nbits + 16
        for i in range(0,maxi,2):
            result = result + term / (i+1)
            term = term*x2*(i+1) / (i+2)
            if abs(term) < test:
                break
        if i >= (maxi-2):
            raise ArithmeticError, "longfloat asin did not converge"

        if didpi6:
            result = result + pi/6
            
        if phase < 0:
            result = -result
        result = extra + result

        longfloat.nbits = longfloat.nbits - 8
        result.truncate()
        if input_phase < 0:
            result = -result
        
        return result

    def atan(x):
        '''
        atan(x) = sum(n) (-1)^n x^(2n+1) / (2n+1)
        %       = x - x^3/3 + x^5/5 - x^7/7      |x| < 1

        atan(x) = - atan(-x) for all x

        atan(x) = pi/2 - atan(1/x) for x >= 0

        atan(x) = pi/2 - acos(x/sqrt(x^2+1)) for all x

        atan(x) = asin(x/sqrt(x^2+1)) for all x

        atan(x) = acot(-x) - pi/2 = pi/2 - acot(x) for all x

        atan(x) + atan(y) = atan((x+y)/(1-x*y))      
        '''

    def atan2(self,x,y):
        '''
        atan(y/x) with apropriate phase
        '''

    def eps(self):
        '''
        Returns eps the smallest number such that 1+eps > 1.
        '''
        eps = longfloat((-longfloat.nbits,1))

        #print ' eps   ', repr(eps)
        #print ' eps   ', eps
        #print ' 1+eps ', 1+eps
        #print ' 1+eps/2', 1+longfloat((-longfloat.nbits-1,1))
        #print ' 1+eps/4', 1+longfloat((-longfloat.nbits-2,1))

        return eps
        
    def __nbits_in_mantissa(self):
        ''' Return the no. of bits in the mantissa '''
        
        absm = abs(self.mantissa)
        if absm == 0:
            return 0

        # Find bounds
        klo = khi = 2*self.nbits+2
        if (absm >> khi):
            khi = khi + 1
            while (absm >> khi):
                klo = khi
                khi = khi << 1
        else:
            klo = klo - 2
            while not (absm >> klo):
                khi = klo
                klo = klo >> 1
            
        # Binary search
        delta = (khi - klo) >> 1
        while delta > 0:
            k = klo + delta
            if (absm >> k):
                klo = k
            else:
                khi = k
            delta = (khi - klo) >> 1
        return klo
        
    def truncate(self):
        '''
        Truncates self to the required number of significant binary digits.
        '''
        nchopped = self.__nbits_in_mantissa() - longfloat.nbits
        if nchopped > 0:
            round = (self.mantissa > 0) and ((self.mantissa >> (nchopped-1)) & 1L)
            self.mantissa = self.mantissa >> nchopped
            if round: self.mantissa = self.mantissa + 1
            self.exponent = self.exponent + nchopped
            if not self.mantissa: self.exponent = 0
 
    def other_truncate(self, nbits=None):
        '''
        Truncates self to the required number of significant binary digits.
        '''
        if nbits == None:
            nbits = self.nbits
        
        eold = e = self.exponent
        m = self.mantissa

        k = nbits  # Multiplications typically double the no. of bits
        test  = (1L << nbits)-1  # The biggest permisible mantissa

        flag = (abs(m) > test)
        while flag:
            mk = m >> k
            flag2 = (abs(mk) > test)
            if (k == 1) or flag2:
                m = mk
                flag = flag2
                e = e + k
                if k > 1:
                    k = k >> 1
            else:
                k = k >> 1

        if self.debug:
            print 'longfloat: truncate: nbits lost =', e - eold
        if m == 0L: e = 0

        # Comment out the next three lines to stop rounding
        if (e-eold):
            if (self.mantissa >> (e-eold-1)) & 1:
                m = m + 1

        self.exponent = e
        self.mantissa = m

    def __div__(self, other):
        return self*other.reciprocal()

    def __rdiv__(self,other):
        return other*self.reciprocal()
        

def __mp_do_test(test_nbits,ntest=10,ntime=10):
    from whrandom import random, seed
    from sys import exit
    from time import clock

    seed((3*test_nbits-1)%256,(5*test_nbits-2)%256,(7*test_nbits-3)%256)

    savenbits = longfloat.nbits
    longfloat.nbits = test_nbits
    print '\n\n Testing extended precision routines with nbits =', \
          longfloat.nbits
    print ' '
    
    def nbits(eps):
        return int(-log(abs(eps))/log(2.0))

    fac = longfloat(1)/longfloat(3)
    
    # Test reciprocal
    nbits_worst = longfloat.nbits
    n = ntest
    print ' n =',n
    maxx = 0
    minx = minabsx = 1e300
    while n > 0:
        x = longfloat((random()-0.5)*10.0)
        x = x*x*x*x*x # To fill out the number of bits
        if x:
            maxx = max(x,maxx)
            minx = min(x,minx)
            minabsx = min(abs(x),minabsx)
            eps = abs(float(x/x - 1.0))
            if eps:
                nbits_worst = min(nbits_worst, nbits(eps))
        n = n - 1
    print (' maxx=%+9.1e minx=%+9.1e minabsx=%+9.1e ' % (maxx, minx, minabsx))
    print ' Reciprocal:  eps = (x(1/x)-1)       :  worst nbits = ', nbits_worst

    x = longfloat(33)
    start = clock()
    n = ntime
    while n > 0:
        x.reciprocal()
        n = n - 1
    used = (clock() - start)/ntime
    print ' Time per call ', used
    print ' '
    
    # Test sqrt
    nbits_worst = longfloat.nbits
    n = ntest
    print ' n =',n
    maxx = 0
    minx = minabsx = 1e300
    while n > 0:
        x = longfloat(random()*10.0)
        x = x*x*x*x*x # To fill out the number of bits
        if x:
            maxx = max(x,maxx)
            minx = min(x,minx)
            minabsx = min(abs(x),minabsx)
            y = x.sqrt()
            # If there is a relative error eps in y=sqrt(x) then
            # eps = (y*y - x)/(2x)
            eps = 0.5*abs(float((y*y - x))/float(x))
            if eps:
                nbits_worst = min(nbits_worst, nbits(eps))
        n = n - 1
    print (' maxx=%+9.1e minx=%+9.1e minabsx=%+9.1e ' % (maxx, minx, minabsx))
    print ' Square root: eps = (sqrt(x)^2-x)/2x :  worst nbits = ', nbits_worst
    x = longfloat(33)
    start = clock()
    n = ntime
    while n > 0:
        x.rsqrt()
        n = n - 1
    used = (clock() - start)/ntime
    print ' Time per call ', used
    print ' '
    
    # Test log(exp)
    nbits_worst = longfloat.nbits
    n = ntest
    print ' n =',n
    maxx = 0
    minx = minabsx = 1e300
    while n > 0:
        x = longfloat((random()-0.5)*10.0)
        x = x*x*x*x*x # To fill out the number of bits
        if x:
            maxx = max(x,maxx)
            minx = min(x,minx)
            minabsx = min(abs(x),minabsx)
            y = x.exp()
            eps = abs(float((y.log() - x))/float(x))
            if eps:
                # If there is a relative error eps1 in log(x)
                # and eps2 in exp(x), then the relative error in
                # -- exp(log(x)) is eps2 + eps1*|logx|
                # -- log(exp(x)) is eps1 + eps2/|x|
                absx = abs(float(x))
                eps = eps /(1.0 + 1.0/absx)
                nbits_worst = min(nbits_worst, nbits(eps))
        n = n - 1
    print (' maxx=%+9.1e minx=%+9.1e minabsx=%+9.1e ' % (maxx, minx, minabsx))
    print ' log & exp :  eps = (log(exp(x))-x)/x:  worst nbits = ', nbits_worst
    x = longfloat(33)
    scale = longfloat(77)
    start = clock()
    n = ntime
    while n > 0:
        x.log()
        n = n - 1
        x = x * scale
    used = (clock() - start)/ntime
    print ' Time per call for log ', used
    x = longfloat(3)
    start = clock()
    n = ntime
    while n > 0:
        x.exp()
        n = n - 1
    used = (clock() - start)/ntime
    print ' Time per call for exp ', used
    print ' '

    # Test exp(log)
    nbits_worst = longfloat.nbits
    n = ntest
    print ' n =',n
    maxx = 0
    minx = minabsx = 1e300
    while n > 0:
        x = longfloat(random()*10.0)
        x = x*x*x*x*x # To fill out the number of bits
        if x:
            maxx = max(x,maxx)
            minx = min(x,minx)
            minabsx = min(abs(x),minabsx)
            y = x.log()
            eps = abs(float((y.exp() - x))/float(x))
            if eps:
                abslogx = abs(float(y))
                # If there is a relative error eps1 in log(x)
                # and eps2 in exp(x), then the relative error in
                # -- exp(log(x)) is eps2 + eps1*|logx|
                # -- log(exp(x)) is eps1 + eps2/|x|
                eps = eps / (1.0 + abslogx)
                nbits_worst = min(nbits_worst, nbits(eps))
        n = n - 1
    print (' maxx=%+9.1e minx=%+9.1e minabsx=%+9.1e ' % (maxx, minx, minabsx))
    print ' log & exp :  eps = (exp(log(x))-x)/x:  worst nbits = ', nbits_worst
    print ' '

    # Test sine and asin
    # Note that asin(sin(x)) is inherently unstable so test just sin(asin(x))
    nbits_worst = longfloat.nbits
    n = ntest
    print ' n =',n
    maxx = 0
    minx = minabsx = 1e300
    xworst = 0.0
    while n > 0:
        x = longfloat((random()-0.5))*2.0  # x in [-1,1]
        if x:
            maxx = max(x,maxx)
            minx = min(x,minx)
            minabsx = min(abs(x),minabsx)
            eps = float(x.asin().sin() - x)/float(x)
            if eps:
                nbb = nbits(eps)
                if nbb < nbits_worst:
                    xworst = x
                    nbits_worst = nbb
                #nbits_worst = min(nbits_worst, nbits(eps))
                #print float(x), nbits(eps)
        n = n - 1
    print (' maxx=%+9.1e minx=%+9.1e minabsx=%+9.1e ' % (maxx, minx, minabsx))
    print ' sin       : eps = (sin(asin(x))-x)/x:  worst nbits = ', nbits_worst
    print ' xworst ', xworst
    x = longfloat(2.7)
    start = clock()
    n = ntime
    while n > 0:
        x.sin()
        n = n - 1
    used = (clock() - start)/ntime
    print ' Time per call for sin  ', used
    x = longfloat(0.33333)
    start = clock()
    n = ntime
    while n > 0:
        x.asin()
        n = n - 1
    used = (clock() - start)/ntime
    print ' Time per call for asin ', used

    
    # Test cosine and acos
    nbits_worst = longfloat.nbits
    n = ntest
    print ' n =',n
    maxx = 0
    minx = minabsx = 1e300
    xworst = 0.0
    while n > 0:
        x = longfloat((random()-0.5))*2.0  # x in [-1,1]
        if x:
            maxx = max(x,maxx)
            minx = min(x,minx)
            minabsx = min(abs(x),minabsx)
            eps = float(x.acos().cos() - x)/float(x)
            if eps:
                nbb = nbits(eps)
                if nbb < nbits_worst:
                    xworst = x
                    nbits_worst = nbb
        n = n - 1
    print (' maxx=%+9.1e minx=%+9.1e minabsx=%+9.1e ' % (maxx, minx, minabsx))
    print ' cos       : eps = (cos(acos(x))-x)/x:  worst nbits = ', nbits_worst
    print ' xworst ', xworst
    x = longfloat(2.7)
    start = clock()
    n = ntime
    while n > 0:
        x.cos()
        n = n - 1
    used = (clock() - start)/ntime
    print ' Time per call for cos  ', used
    x = longfloat(0.33333)
    start = clock()
    n = ntime
    while n > 0:
        x.acos()
        n = n - 1
    used = (clock() - start)/ntime
    print ' Time per call for acos ', used

    

    # Test e and pi
    maplee  = '2.718281828459045235360287471352662497757247093699959574966967627724076630353547594571382178525166427427466391932003059921817413596629043572900334295260595630738132328627943490763233829880753195251019011573834187930702154089149934884167509244761460668082264800168477411853742345442437107539077744992069551702761838606261331384583000752044933826560297606737113200709328709127443747047230696977209310141692836819025515108657463772111252389784425056953696770785449969967946864454905987931636889230098793127736178215424999229576351482208269895193668033182528869398496465105820939239829488793320362509443117301238197068416140397019837679320683282376464804295311802328782509819455815301756717361332069811250996181881593041690351598888519345807273866738589422879228499892086805825749279610484198444363463244968487560233624827041978623209002160990235304369941849146314093431738143640546253152096183690888707016768396424378140592714563549061303107208510383750510115747704171898610687396965521267154688957035035'
    maplepi = '3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647093844609550582231725359408128481117450284102701938521105559644622948954930381964428810975665933446128475648233786783165271201909145648566923460348610454326648213393607260249141273724587006606315588174881520920962829254091715364367892590360011330530548820466521384146951941511609433057270365759591953092186117381932611793105118548074462379962749567351885752724891227938183011949129833673362440656643086021394946395224737190702179860943702770539217176293176752384674818467669405132000568127145263560827785771342757789609173637178721468440901224953430146549585371050792279689258923542019956112129021960864034418159813629774771309960518707211349999998372978049951059731732816096318595024459455346908302642522308253344685035261931188171010003137838752886587533208381420617177669147303598253490428755468731159562863882353787593751957781857780532171226806613001927876611195909216420199'
    ndigit = int(longfloat.nbits/(log(10)/log(2)))
    ndigit = min(1001,ndigit+2)
    
    print (' The error in e and pi should be in the last d.p. (about %9.1e)' %
           (3.0*(0.5**longfloat.nbits)))
    print '    Computed e   ', longfloat(1).exp()
    print '       Maple e    ', maplee[:ndigit]
    print '    Computed pi  ', longfloat(1).pi()
    print '       Maple pi   ', maplepi[:ndigit]
    print ' ' 
    print ' nbits at end ', longfloat.nbits
    print ' '
    
    longfloat.nbits = savenbits
    print ' restored longfloat.nbits ', savenbits

if __name__ == "__main__":
    from time import clock
    I = interpolators()
##     print ' EXP '
##     vv = I.exp(longfloat(0.25))
##     start = clock()
##     for i in range(10):
##         vv = I.exp(longfloat(0.25))
##     used = (clock() - start)/10
##     print ' NEW ', used
##     start = clock()
##     for i in range(10):
##         vv = longfloat(0.25).exp()
##     used = (clock() - start)/10
##     print ' OLD ', used
    ## print ' SIN '
##     print I.sin(0.25,1e-16)
##     print sin(0.25)
##     print ' COS '
##     print I.cos(0.25,1e-16)
##     print cos(0.25)
    print ' LOG '
    print I.log(longfloat(0.9))
    print longfloat(0.9).log()
    q = longfloat(0.72)
    I.log(q)
##     start = clock()
##     for i in range(10):
##         vv = I.log(q)
##     used = (clock() - start)/10
##     print ' NEW ', used
##     start = clock()
##     for i in range(10):
##         vv = q.log()
##     used = (clock() - start)/10
##     print ' OLD ', used
    import profile
    profile.run('for i in range(100): q.reciprocal()')
##     print ' ASIN by inversion'
##     print I.invert(0.25,I.sindsin,0.3,1e-16)
##     print asin(0.25)
##     print ' ACOS by inversion'
##     print I.invert(0.25,I.cosdcos,1.1,1e-16)
##     print acos(0.25)
##     print ' EXP=ALOG by inversion'
##     print I.invert(0.25,I.logdlog,1.0,1e-16)
##     print exp(0.25)
##     print ' LOG=AEXP by inversion'
##     print I.invert(1.35,I.expdexp,0.2,1e-16)
##     print log(1.35)
    #print ' HI'
    #__mp_do_test(104,ntest=3,ntime=3)
    #__mp_do_test(208,ntest=3,ntime=3)
    #__mp_do_test(500)
    #__mp_do_test(1000)
else:
    #print ' reimported mp'
    pass

## def sqrt(x):
##     return longfloat(x).sqrt()

## def sin(x):
##     return longfloat(x).sin()

## def cos(x):
##     return longfloat(x).cos()

## def asin(x):
##     return longfloat(x).asin()

## def acos(x):
##     return longfloat(x).acos()

## def exp(x):
##     return longfloat(x).exp()

## def log(x):
##     return longfloat(x).log()

