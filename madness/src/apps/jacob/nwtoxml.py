import sys

basis_name = "6-31G"

print '<?xml version="1.0" ?>'
print '<name>'
print "  ", basis_name
print '</name>'

while 1:
    line = sys.stdin.readline()
    if not line: break

    if line[0:5] == "basis":
        symbol = line.split()[1].split("_")[0][1:]
        print '<basis symbol="%s">' % symbol

        line = sys.stdin.readline().split()
        while 1:
            if line[0] == "end":
                break
            shelltype = line[1]
            coeff = []
            pcoeff = []
            expnt = []
            while 1:
                line = sys.stdin.readline().split()
                e,c,p = 0.0,0.0,0.0
                try:
                    e,c = float(line[0]), float(line[1])
                    if shelltype == 'SP': p = float(line[2])
                    expnt.append(e)
                    coeff.append(c)
                    pcoeff.append(p)
                except:
                    break
            nprim = len(expnt)
            if shelltype == "SP": shelltype = "L"
            print '  <shell type="%s" nprim="%d">' % (shelltype, nprim)
            print '    <exponents>'
            for e in expnt: print "      %.8f" % e
            print '    </exponents>'
            if shelltype == 'L':
                print '    <scoefficients>'
                for e in coeff: print "      %.8f" % e
                print '    </scoefficients>'
                print '    <pcoefficients>'
                for e in pcoeff: print "      %.8f" % e
                print '    </pcoefficients>'
            else:
                print '    <coefficients>'
                for e in coeff: print "      %.8f" % e
                print '    </coefficients>'

            print '  </shell>'

        print '</basis>'

