#
#@BEGIN LICENSE
#
# PSI4: an ab initio quantum chemistry software package
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
#@END LICENSE
#

from __future__ import print_function
"""Queries the PubChem database using a compound name (i.e. 1,3,5-hexatriene)
   to obtain a molecule string that can be passed to Molecule. ::

      results = getPubChemObj("1,3,5-hexatriene")

      Results is an array of results from PubChem matches to your query.
        for entry in results:
           entry["CID"]         => PubChem compound identifer
           entry["IUPAC"]       => IUPAC name for the resulting compound
           entry["PubChemObj"]  => instance of PubChemObj for this compound

           entry["PubChemObj"].getMoleculeString()   => returns a string compatible
                                                        with PSI4's Molecule creation

"""
try:
    # Python 2 syntax
    from urllib2 import urlopen
    from urllib2 import quote
    from urllib2 import URLError
except ImportError:
    # Python 3 syntax
    from urllib.request import urlopen
    from urllib.parse import quote
    from urllib.error import URLError
import xml.etree.ElementTree as ET
import time
import gzip
import re
import sys
import os


class PubChemObj(object):

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
        """Function to return the SDF (structure-data file) of the PubChem object."""
        if (len(self.dataSDF) == 0):
            def extract_xml_keyval(xml, key):
                """ A useful helper function for parsing a single key from XML. """
                matches = list(xml.iter(key))
                if len(matches) == 0:
                    return None
                elif len(matches) == 1:
                    return matches[0].text
                else:
                    raise TooManyMatchesException

            url = "http://pubchem.ncbi.nlm.nih.gov/pug/pug.cgi"
            initial_request = """
            <PCT-Data>
              <PCT-Data_input>
                <PCT-InputData>
                  <PCT-InputData_download>
                    <PCT-Download>
                      <PCT-Download_uids>
                        <PCT-QueryUids>
                          <PCT-QueryUids_ids>
                            <PCT-ID-List>
                              <PCT-ID-List_db>pccompound</PCT-ID-List_db>
                              <PCT-ID-List_uids>
                                <PCT-ID-List_uids_E>%d</PCT-ID-List_uids_E>
                              </PCT-ID-List_uids>
                            </PCT-ID-List>
                          </PCT-QueryUids_ids>
                        </PCT-QueryUids>
                      </PCT-Download_uids>
                      <PCT-Download_format value="sdf"/>
                    </PCT-Download>
                  </PCT-InputData_download>
                </PCT-InputData>
              </PCT-Data_input>
            </PCT-Data>
            """ % self.cid

            print("\tPosting PubChem request for CID %d" % self.cid)
            server_response = urlopen(url, initial_request).read()
            xml = ET.fromstring(server_response)
            for attempt in range(4):
                if attempt == 3:
                    raise PubChemTimeoutError
                # Look for a download location in the XML response
                download_url = extract_xml_keyval(xml, 'PCT-Download-URL_url')
                if download_url:
                    print("\tDownloading from PubChem...")
                    # We have a download location; gogetit
                    file_name = '_psi4_pubchem_temp_download_.tgz'
                    response = urlopen(download_url)
                    data = response.read()
                    fh = open(file_name, 'wb')
                    fh.write(data)
                    fh.close()
                    unzipped = gzip.open(file_name)
                    self.dataSDF = unzipped.read()
                    unzipped.close()
                    try:
                        os.remove(file_name)
                    except:
                        pass
                    print("\tDone!")
                    break
                # We didn't find a download location yet.
                if attempt == 0:
                    # If this is the first try, take a ticket
                    ticket = extract_xml_keyval(xml, 'PCT-Waiting_reqid')
                    #print("ticket = " + ticket)
                    statusrequest = """<PCT-Data>
                      <PCT-Data_input>
                        <PCT-InputData>
                          <PCT-InputData_request>
                            <PCT-Request>
                              <PCT-Request_reqid>%d</PCT-Request_reqid>
                              <PCT-Request_type value="status"/>
                            </PCT-Request>
                          </PCT-InputData_request>
                        </PCT-InputData>
                      </PCT-Data_input>
                    </PCT-Data>
                    """ % int(ticket)
                if ticket:
                    # Wait 10 seconds...
                    print("\tPubChem result not available yet, will try again in 10 seconds...")
                    time.sleep(10)
                    # ...and ask for an update on the progress
                    server_response = urlopen(url, statusrequest).read()
                    xml = ET.fromstring(server_response)
                else:
                    # We can't find a ticket number, or a download location. Bail.
                    raise PubChemDownloadError
        return self.dataSDF

    def name(self):
        """Function to return the IUPAC name of the PubChem object."""
        return self.iupac

    def getCartesian(self):
        """Function to return a string of the atom symbol and XYZ
        coordinates of the PubChem object.

        """
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
        """Function to obtain preferentially a molecule string
        through getCartesian() or a query string otherwise.

        """
        try:
            temp = self.getCartesian()
        except Exception as e:
            raise
        molstr = "%d\n%s\n%s" % (self.natom, self.iupac, temp)
        return molstr

    def getMoleculeString(self):
        """Function to obtain a molecule string through
        getCartesian() or fail.
        """
        try:
            return self.getCartesian()
        except Exception as e:
            return e.message


def getPubChemResults(name):
    """Function to query the PubChem database for molecules matching the
    input string. Builds a PubChem object if found.

    """
    url = 'http://www.ncbi.nlm.nih.gov/sites/entrez?db=pccompound&term=%s&format=text' % quote(name)
    print("\tSearching PubChem database for %s" % (name))
    try:
        loc = urlopen(url)
    except URLError as e:
        msg = "\tPubchemError\n%s\n\treceived when trying to open\n\t%s\n" % (str(e), url)
        msg += "\tCheck your internet connection, and the above URL, and try again.\n"
        raise Exception(msg)
    data = loc.read()

    ans = []
    l = data.find(b"<pre>")
    l = data.find(b"\n", l)
    i = 1
    while(True):
        l = data.find(str("%d. " % i).encode(sys.getdefaultencoding()), l)
        if l == -1:
            break
        tag = b"MF: "
        l = data.find(tag, l) + len(tag)
        mf = data[l:data.find(b'\n', l)].decode(sys.getdefaultencoding())
        tag = b"IUPAC name: "
        l = data.find(tag, l) + len(tag)
        iupac = data[l:data.find(b'\n', l)].decode(sys.getdefaultencoding())
        tag = b"CID:"
        l = data.find(tag, l) + len(tag)
        #if l == 4:
        #    break
        cid = int(data[l:data.find(b"\n", l)])
        l = data.find(b'\t', l) + 1

        pubobj = PubChemObj(cid, mf, iupac)
        ans.append(pubobj)
        i += 1

    print("\tFound %d results" % (len(ans)))
    return ans

if __name__ == "__main__":
    try:
        obj = getPubChemResults("1-methoxy-4-[(E)-prop-1-enyl]benzene")
        #obj = getPubChemResults("sodium benzenesulfonate")
    except Exception as e:
        print(e.message)

    for r in obj:
        print(r)
        print(r.getMoleculeString())
