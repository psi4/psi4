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

import urllib2
import re

class PubChemObj:

    def __init__(self, cid, mf, iupac):
        self.url = 'http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi'
        self.cid = cid
        self.mf = mf
        self.iupac = iupac
        self.natom = 0
        self.dataSDF = ''

    def __str__(self):
        return "%17d   %s\n" % (self.cid, self.iupac)

    def getSDF(self):
        if (len(self.dataSDF) == 0):
            # When completed uncomment the following:
            url = self.url+'?cid='+urllib2.quote(str(self.cid))+'&disopt=3DDisplaySDF'
            try:
                location = urllib2.urlopen(url)
            except urllib2.URLError, e:
                msg =  "\tPubchemError\n%s\n\treceived when trying to open\n\t%s\n" % (str(e),  url)
                msg += "\tCheck your internet connection, and the above URL, and try again.\n"
                raise Exception(msg)
            self.dataSDF = location.read()
            #f = open("TEST", "w")
            #f.write(self.dataSDF)
        return self.dataSDF

    def name(self):
        return self.iupac

    def getCartesian(self):
        try:
            sdfText = self.getSDF()
        except Exception as e:
            raise e

        # Find
        # NA NB                        CONSTANT
        # 14 13  0     0  0  0  0  0  0999 V2000
        m = re.search(r'^\s*(\d+)\s+(?:\d+\s+){8}V2000$', sdfText, re.MULTILINE)
        self.natom = 0
        if (m):
            self.natom = int(m.group(1))

        if (self.natom == 0):
            raise Exception("PubchemError\n Cannot find the number of atoms.  3D data doesn't appear\n" +
                            "to be available for %s.\n" % self.iupac)

        lines = re.split('\n', sdfText)

        #  3.7320   -0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        NUMBER = "((?:[-+]?\\d*\\.\\d+(?:[DdEe][-+]?\\d+)?)|(?:[-+]?\\d+\\.\\d*(?:[DdEe][-+]?\\d+)?))"
        atom_re = re.compile(r'^\s*' + NUMBER + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*(\w+)(?:\s+\d+){12}')

        molecule_string = "PubchemInput\n"

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
        try:
            temp = self.getCartesian()
        except Exception as e:
            raise
        molstr = "%d\n%s\n%s" % (self.natom, self.iupac, temp)
        return molstr

    def getMoleculeString(self):
        try:
            return self.getCartesian()
        except Exception as e:
            return e.message
   
def getPubChemResults(name):
    url = 'http://www.ncbi.nlm.nih.gov/sites/entrez?db=pccompound&term=%s&format=text' % (urllib2.quote(name))
    print "\tSearching PubChem database for %s" % (name)
    try:
        loc = urllib2.urlopen(url)
    except urllib2.URLError as e:
        msg =  "\tPubchemError\n%s\n\treceived when trying to open\n\t%s\n" % (str(e),  url)
        msg += "\tCheck your internet connection, and the above URL, and try again.\n"
        raise Exception(msg)
    data = loc.read()

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

    print "\tFound %d results" % (len(ans))
    return ans

if __name__ == "__main__":
    try:
        obj = getPubChemResults("1-methoxy-4-[(E)-prop-1-enyl]benzene")
        #obj = getPubChemResults("sodium benzenesulfonate")
    except Exception as e:
        print e.message

    for r in obj:
        print r
        print r.getMoleculeString()

