import sys
from array import array

__kmax = 60
__asciifile = "coeffs"
__checksum = (1.57015282218922890,1.48602374574943450)

def matrix(n,m): 
    a = range(n) 
    for i in range(n): 
        a[i] = array('d',[0.0]*m) 
    return a 

def writemat(k,file,a):
    for i in range(k):
        for j in range(k):
            if a[i][j] == 0.0:
                file.write("0.0\n")
            else:
                file.write("%.17e\n" % a[i][j])

def readmat(k,file):
    a = matrix(k,k)
    for i in range(k):
        for j in range(k):
            a[i][j] = float(file.readline())

    return a

def writeascii(kmax,filename):
    file = open(filename,'w')
    for k in range(1,kmax+1):
        h0,h1,g0,g1 = twoscalecoeffs(k)
        writemat(k,file,h0)
        writemat(k,file,g0)

def readascii(kmax,filename):
    global __cache
    __cache = {}
    file = open(filename,'r')
    for k in range(1,kmax+1):
        h0 = readmat(k,file)
        g0 = readmat(k,file)
        h1 = matrix(k,k)
        g1 = matrix(k,k)
        for i in range(k):
            for j in range(k):
                hphase = (-1)**(i+j)
                gphase = (-1)**(i+j+k)
                h1[i][j] = h0[i][j]*hphase
                g1[i][j] = g0[i][j]*gphase

        __cache[k] = h0,h1,g0,g1

    sum0, sum1 = checksum(kmax)
    if abs(__checksum[0]-sum0) > 1e-16 or \
       abs(__checksum[1]-sum1) > 1e-16:
        print "%.17e %.17e" % __checksum
        print "%.17e %.17e" % (sum0, sum1)
        raise ValueError,"two-scale checksum?"
    
def checksum(kmax):
    sum0 = sum1 = 0.0
    for k in range(1,kmax+1):
        h0,h1,g0,g1 = __cache[k]
        for i in range(k):
            for j in range(k):
                ij = float(j+k)/(i+k)
                sum0 = sum0 + ij*(abs(h0[i][j]) + abs(h1[i][j]) +
                                  abs(g0[i][j]) + abs(g1[i][j]))
                sum1 = sum1 + ij*(h0[i][j] + h1[i][j] +
                                  g0[i][j] + g1[i][j])
                while sum0 > 2.0: sum0 = sum0 * 0.25
                while abs(sum1) > 2.0: sum1 = sum1 * 0.25

    return (sum0, sum1)

if __name__ == "__main__":
    from twoscalecoeffs import twoscalecoeffs
    writeascii(__kmax,__asciifile)
    readascii(__kmax,__asciifile)
