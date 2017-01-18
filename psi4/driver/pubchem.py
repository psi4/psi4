#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2017 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
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
# @END LICENSE
#

from __future__ import print_function,absolute_import
"""Queries the PubChem database using a compound name (i.e. 1,3,5-hexatriene)
   to obtain a molecule string that can be passed to Molecule. ::

      results = getPubChemObj("1,3,5-hexatriene")

      Results is an array of results from PubChem matches to your query.
        for entry in results:
           entry["CID"]         => PubChem compound identifer
           entry["IUPAC"]       => IUPAC name for the resulting compound
           entry["PubChemObj"]  => instance of PubChemObj for this compound

           entry["PubChemObj"].getMoleculeString()   => returns a string compatible
                                                        with Psi4's Molecule creation

"""

try:
    # Python 2 syntax
    from urllib2 import urlopen, Request
    from urllib2 import quote
    from urllib2 import URLError
except ImportError:
    # Python 3 syntax
    from urllib.request import urlopen, Request
    from urllib.parse import quote
    from urllib.error import URLError
import xml.etree.ElementTree as ET
import json
import time
import gzip
import re
import sys
import os
from psi4.driver.p4util.exceptions import *

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
                try:
                    # Python 2.6 (ElementTree 1.2 API)
                    matches = xml.getiterator(key)
                except:
                    # Python 2.7 (ElementTree 1.3 API)
                    matches = list(xml.iter(key))
                if len(matches) == 0:
                    return None
                elif len(matches) == 1:
                    return matches[0].text
                else:
                    print(matches)
                    raise ValidationError("""PubChem: too many matches found %d""" % (len(matches)))

            url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%d/SDF?record_type=3d' % self.cid
            req = Request(url, headers={'Accept' : 'chemical/x-mdl-sdfile'})
            try:
                self.dataSDF = urlopen(req).read().decode('utf-8')
            except URLError as e:
                msg = "Unable to open\n\n%s\n\ndue to the error\n\n%s\n\n" %(url, e)
                msg += "It is possible that 3D information does not exist for this molecule in the PubChem database\n"
                print(msg)
                raise ValidationError(msg)
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

        if self.natom == 0:
            raise ValidationError("PubChem: Cannot find the number of atoms.  3D data doesn't appear\n" +
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
    print("\tSearching PubChem database for %s" % (name))
    url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/%s/property/IUPACName,MolecularFormula/JSON' % quote(name)
    try:
        response = urlopen(url)
    except URLError as e:
        msg = "\tPubchemError\n%s\n\treceived when trying to open\n\t%s\n" % (str(e), url)
        msg += "\tCheck your internet connection, and the above URL, and try again.\n"
        raise ValidationError(msg)
    data = json.loads(response.read().decode('utf-8'))
    results = []
    for d in data['PropertyTable']['Properties']:
        pubobj = PubChemObj(d['CID'], d['IUPACName'], d['IUPACName'])
        results.append(pubobj)

    print("\tFound %d result%s" % (len(results), "" if len(results)==1 else "s"))
    return results


if __name__ == "__main__":
    try:
        #obj = getPubChemResults("1-methoxy-4-[(E)-prop-1-enyl]benzene")
        #obj = getPubChemResults("sodium benzenesulfonate")
        obj = getPubChemResults("4-[bis(4-hydroxyphenyl)methyl]phenol")
    except Exception as e:
        print(e.message)

    for r in obj:
        print(r)
        print(r.getMoleculeString())
