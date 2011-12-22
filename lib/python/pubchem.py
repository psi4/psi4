# -*- coding: utf-8 -*-

#    OriginDB
#    Copyright (C) 2009  Jason Power <power.jg@gmail.com>
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License along
#    with this program; if not, write to the Free Software Foundation, Inc.,
#    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

import urllib
import pybel

class PubChemObj:

        def __init__(self, cid=0):

                self.url = 'http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi'
                self.cid = cid
                self.dataXML = ''
                self.dataSDF = ''


        def getXML(self):

                location = urllib.urlretrieve(self.url+'?cid='+str(self.cid)+'&disopt=SaveXML')[0]

                dataFile = open(location, 'r')
                self.dataXML = dataFile.read()

        def writeXML(self, filename):

                if not self.dataXML:
                        self.getXML()
                if not self.dataXML:
                        return

                f = open(filename, 'w')
                f.write(self.dataXML)


        def getSDF(self):

                location = urllib.urlretrieve(self.url+'?cid='+str(self.cid)+'&disopt=SaveSDF')[0]

                dataFile = open(location, 'r')
                self.dataSDF = dataFile.read()

                return self.dataSDF

        def getMol(self):

                return pybel.readstring('sdf', self.dataSDF)


def searchPubChem(name):
        url = "http://www.ncbi.nlm.nih.gov/sites/entrez?db=pccompound&term=%s&format=text" % (name)

        loc = urllib.urlretrieve(url)[0]
        dataFile = open(loc, 'r')
        data = dataFile.read()

        ans = {}

        l = data.find("<pre>")
        l = data.find("\n", l)
        for i in range(1, 21):
                l = data.find("%s:" % (i), l)
                l = data.find("CID: ", l) + 5
                if l == 4:
                        break
                cid = int(data[l:data.find("\n", l)])
                l = data.find("\t", l)+1
                name = data[l:data.find(";", l)]
                if name.find("\n") != -1:
                        name = data[l:data.find("\n", l)]

                ans[cid] = name

        return ans

