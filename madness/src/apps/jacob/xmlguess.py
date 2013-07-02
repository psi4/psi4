import sys, string

data = sys.stdin
basis = "STO-3G"

def find_string(s):
    n = len(s)
    while 1:
        line = data.readline()
        if not line:
            sys.exit(0)

        if line[0:n] == s:
            return line

def read_geometry():
    find_string("            XYZ format geometry")

    data.readline()
    data.readline()
    data.readline()
    atom = data.readline().split()[0]
    return atom

def read_numao():
    line = find_string("          AO basis - number of functions:").split()
    n = int(line[6])
    return n

def read_energy():
    find_string("     <<<<<<              called scfcvg")
    line = find_string("         Total DFT energy =").split()
    energy = float(line[4])
    return energy

def read_occn(n, ab):
    occ = [0]*n
    eps = [0]*n
    test = "                    DFT Final %s Molecular Orbital Analysis" % ab
    if ab == "Beta":
        test = " " + test
    find_string(test)
    
    for i in range(n):
        line = find_string(" Vector ").split()
        occ[i] = float(string.replace(line[2].split("=")[1],"D","e"))
        e = line[4]
        if e[0] == "S":
            e = line[3].split("=")[1]
        eps[i] = float(string.replace(e,"D","e"))
        
    return (occ,eps)

def matrix(n,m):
    V = [0]*n
    for i in range(n):
        V[i] = [0]*m
    return V

def read_matrix(n,m):
    V = matrix(n,m)
    for jlo in range(0,m,6):
        jhi = min(m,jlo+6)
        data.readline() # blank
        data.readline() # numbers
        data.readline() # underscore
        for i in range(n):
            line = data.readline().split()
            for jj in range(jhi-jlo):
                j = jlo + jj
                V[i][j] = float(line[1+jj])
    return V

def read_orbitals(n):
    find_string(" global array:")
    V = read_matrix(n,n)
    return V

def make_density(n, aocc, amo, bocc, bmo):
    D = matrix(n,n)
    for i in range(n):
        for j in range(n):
            sum = 0.0
            for k in range(n):
                sum = sum + aocc[k]*amo[i][k]*amo[j][k] + bocc[k]*bmo[i][k]*bmo[j][k]
            D[i][j] = sum
    return D

def print_matrix(n, m, a):
    for i in range(n):
        for j in range(m):
            print "%12.5f" % a[i][j],
        print ""

def print_vector(n, a):
    for i in range(n):
        print "%12.5f" % a[i],
    print ""

while 1:
    atom = read_geometry()
    n =  read_numao()
    e = read_energy()

    aocc, aeps = read_occn(n, 'Alpha')
    bocc, beps = read_occn(n, 'Beta')

    amo = read_orbitals(n)
    bmo = read_orbitals(n)

    density = make_density(n, aocc, amo, bocc, bmo)

    print "<atomicguess symbol=\"%s\" basis=\"%s\">" % (atom,basis)
    print "   <guessdensitymatrix>"
    print_matrix(n, n, density)
    print "   </guessdensitymatrix>"
    print "   <alphaocc>"
    print_vector(n, aocc)
    print "   </alphaocc>"
    print "   <betaocc>"
    print_vector(n, bocc)
    print "   </betaocc>"
    print "   <alphaeps>"
    print_vector(n, aeps)
    print "   </alphaeps>"
    print "   <betaeps>"
    print_vector(n, beps)
    print "   </betaeps>"
    print "   <alphavectors>"
    print_matrix(n, n, amo)
    print "   </alphavectors>"
    print "   <betavectors>"
    print_matrix(n, n, bmo)
    print "   </betavectors>"
    print "</atomicguess>"
    



    

