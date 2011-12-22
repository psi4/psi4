# PubChem.py
#
# Queries the PubChem database using a compound name (i.e. 1,3,5-hexatriene)
# to obtain a molecule string that can be passed to Molecule.
#
# Example:
#
# results = getPubChemObj("1,3,5-hexatriene")
#
# Results is an array of results from PubChem matches to your query.
#   for entry in results:
#      entry["CID"]         => PubChem compound identifer
#      entry["IUPAC"]       => IUPAC name for the resulting compound
#      entry["PubChemObj"]  => instance of PubChemObj for this compound
#
#      entry["PubChemObj"].getMoleculeString()   => returns a string compatible
#                                                   with PSI4's Molecule creation
#

import urllib
import re

class PubChemObj:

    def __init__(self, cid, mf, iupac):
        self.url = 'http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi'
        self.cid = cid
        self.mf = mf
        self.iupac = iupac
        self.natom = 0
        self.dataSDF = ''

    def getSDF(self):
        if (len(self.dataSDF) == 0):
            print "Querying PubChem for compound %d" % (self.cid)
            # When completed uncomment the following:
            location = urllib.urlretrieve(self.url+'?cid='+str(self.cid)+'&disopt=SaveSDF')[0]
            dataFile = open(location, 'r')
            self.dataSDF = dataFile.read()
            # and remove the following
            #self.dataSDF = sdf

        return self.dataSDF

    def getCartesian(self):
        sdfText = self.getSDF()

        # Find
        # NA NB                        CONSTANT
        # 14 13  0     0  0  0  0  0  0999 V2000
        m = re.search(r'^\s*(\d+)\s+(?:\d+\s+){8}V2000$', sdfText, re.MULTILINE)
        self.natom = 0
        if (m):
            self.natom = int(m.group(1))

        if (self.natom == 0):
            print "Unable to find the number of atoms."
            return None

        lines = re.split('\n', sdfText)

        #  3.7320   -0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        NUMBER = "((?:[-+]?\\d*\\.\\d+(?:[DdEe][-+]?\\d+)?)|(?:[-+]?\\d+\\.\\d*(?:[DdEe][-+]?\\d+)?))"
        atom_re = re.compile(r'^\s*' + NUMBER + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*(\w+)(?:\s+\d+){12}')

        molecule_string = ""

        atom_count = 0
        for line in lines:
            if (not line or line.isspace()):
                continue

            atom_match = atom_re.match(line)
            if atom_match:
                x = float(atom_match.group(1))
                y = float(atom_match.group(2))
                z = float(atom_match.group(3))
                sym = atom_match.group(4)

                atom_count = atom_count + 1

                molecule_string += "%s %10.6f %10.6f %10.6f\n" % (sym, x, y, z)

                if (atom_count == self.natom):
                    break

        return molecule_string

    def getXYZFile(self):
        temp = self.getCartesian()

        molstr = "%d\n%s\n%s" % (self.natom, self.iupac, temp)
        return molstr

    def getMoleculeString(self):
        return self.getCartesian()

def getPubChemResults(name):
    url = 'http://www.ncbi.nlm.nih.gov/sites/entrez?db=pccompound&term=%s&format=text' % (name)

    # When PubChemObj is finished uncomment the following 3 lines
    print "Searching PubChem database for %s" % (name)
    loc = urllib.urlretrieve(url)[0]
    dataFile = open(loc, 'r')
    data = dataFile.read()
    # and remove the following:
    #data = searchresults

    ans = []
    l = data.find("<pre>")
    l = data.find("\n", l)
    for i in range(1, 21):
        l = data.find("%s. " % (i), l)
        if l == -1:
            break
        l = data.find("MF: ", l) + 4
        mf = data[l:data.find("\n", l)]
        l = data.find("IUPAC: ", l) + 7
        iupac = data[l:data.find("\n", l)]
        l = data.find("CID: ", l) + 5
        #if l == 4:
        #    break
        cid = int(data[l:data.find("\n", l)])
        l = data.find("\t", l)+1

        pubobj = PubChemObj(cid, mf, iupac)
        ans.append(pubobj)

    print "Found %d results" % (len(ans))
    return ans

if __name__ == "__main__":
    obj1 = getPubChemResults("1,3,5-hexatriene")
    obj2 = getPubChemResults("benzene")

    print "Results for 1,3,5-hexatriene"
    for mol in obj1:
        print mol["PubChemObj"].getMoleculeString()

    print "Results for benzene"
    for mol in obj2:
        print "IUPAC name: %s\nCompound: %d" % (mol["IUPAC"], mol["CID"])
        print mol["PubChemObj"].getMoleculeString()
