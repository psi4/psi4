# -*- coding: utf-8 -*-
"""
CIRpy

Python interface for the Chemical Identifier Resolver (CIR) by the CADD Group at the NCI/NIH.
https://github.com/mcs07/CIRpy

Originally written by Matt Swain and described in a blog post here:
http://matt-swain.com/post/19633070138/cirpy-a-python-interface-for-the-chemical-identifier

At the point at which I copied this, there was no accompanying license, so I assumed that we
could just copy the code here.  It would not take long to write one for ourselves, but this
appears to be working fine.

Note that some changes have been made to make it (more) compatible with Grendel.  This should
be access through the interface in the web_getter module.
"""


import os
import urllib2
from xml.etree import ElementTree as ET


API_BASE = 'http://cactus.nci.nih.gov/chemical/structure'

available_properties = []

def resolve(input, representation, resolvers=None):
    """ Resolve input to the specified output representation """
    resultdict = query(input, representation, resolvers)
    result = resultdict[0]['value'] if resultdict else None
    if result and len(result) == 1:
        result = result[0]
    return result


def query(input, representation, resolvers=None):
    """ Get all results for resolving input to the specified output representation """
    apiurl = API_BASE+'/%s/%s/xml' % (urllib2.quote(input), representation)
    if resolvers is not None:
        r = ",".join(resolvers)
        apiurl += '?resolver=%s' % r
    result = []
    tree = ET.parse(urllib2.urlopen(apiurl))
    for data in tree.findall(".//data"):
        datadict = {'resolver':data.attrib['resolver'],
                    'notation':data.attrib['notation'],
                    'value':[]}
        for item in data.findall("item"):
            datadict['value'].append(item.text)
        if len(datadict['value']) == 1:
            datadict['value'] = datadict['value'][0]
        result.append(datadict)
    return result if result else None

def get_in_format(input, format='sdf'):
    """ Resolve and get the structure as a str"""
    url = API_BASE+'/%s/file?format=%s' % (urllib2.quote(input), format)
    servefile = urllib2.urlopen(url)
    return str(servefile.read())

def download(input, filename, format='sdf', overwrite=False):
    """ Resolve and download structure as a file """
    url = API_BASE+'/%s/file?format=%s' % (urllib2.quote(input), format)
    servefile = urllib2.urlopen(url)
    if not overwrite and os.path.isfile(filename):
        raise IOError("%s already exists. Use 'overwrite=True' to overwrite it." % filename)
    file = open(filename, "w")
    file.write(servefile.read())
    file.close()


class CIRCachedProperty(object):
    """ Descriptor for caching CIRMolecule properties. """

    def __init__(self, func):
        self._func = func
        self.__name__ = func.__name__
        self.__doc__ = func.__doc__
        globals()["available_properties"].append(func.__name__)

    def __get__(self, obj, obj_class=None):
        if obj is None: return None
        result = obj.__dict__[self.__name__] = self._func(obj)
        return result

class CIRCachedFileFormat(object):
    """ Descriptor for caching CIRMolecule properties. """

    def __init__(self, func):
        self._func = func
        self.__name__ = func.__name__
        self.__doc__ = func.__doc__
        globals()["available_properties"].append(func.__name__)

    def __get__(self, obj, obj_class=None):
        if obj is None: return None
        result = obj.__dict__[self.__name__] = self._func(obj, self.__name__)
        return result


class CIRMolecule(object):
    """Class to hold and cache the structure information for a given CIR input"""

    def __init__(self, input, resolvers=None):
        """ Initialize with a query input """
        self.input = input
        self.resolvers = resolvers

    def __repr__(self):
        return "CIRMolecule(%r, %r)" % (self.input, self.resolvers)


    @CIRCachedProperty
    def stdinchi(self): return resolve(self.input, 'stdinchi', self.resolvers)

    @CIRCachedProperty
    def stdinchikey(self): return resolve(self.input, 'stdinchikey', self.resolvers)

    @CIRCachedProperty
    def smiles(self): return resolve(self.input, 'smiles', self.resolvers)

    @CIRCachedProperty
    def ficts(self): return resolve(self.input, 'ficts', self.resolvers)

    @CIRCachedProperty
    def ficus(self): return resolve(self.input, 'ficus', self.resolvers)

    @CIRCachedProperty
    def uuuuu(self): return resolve(self.input, 'uuuuu', self.resolvers)

    @CIRCachedProperty
    def hashisy(self): return resolve(self.input, 'hashisy', self.resolvers)

    @CIRCachedProperty
    def sdf(self): return resolve(self.input, 'sdf', self.resolvers)

    @CIRCachedProperty
    def names(self): return resolve(self.input, 'names', self.resolvers)

    @CIRCachedProperty
    def iupac_name(self): return resolve(self.input, 'iupac_name', self.resolvers)

    @CIRCachedProperty
    def cas(self): return resolve(self.input, 'cas', self.resolvers)

    @CIRCachedProperty
    def chemspider_id(self): return resolve(self.input, 'chemspider_id', self.resolvers)

    @CIRCachedProperty
    def mw(self): return resolve(self.input, 'mw', self.resolvers)

    @CIRCachedProperty
    def formula(self): return resolve(self.input, 'formula', self.resolvers)

    @CIRCachedProperty
    def h_bond_donor_count(self): return resolve(self.input, 'h_bond_donor_count', self.resolvers)

    @CIRCachedProperty
    def h_bond_acceptor_count(self): return resolve(self.input, 'h_bond_acceptor_count', self.resolvers)

    @CIRCachedProperty
    def h_bond_center_count(self): return resolve(self.input, 'h_bond_center_count', self.resolvers)

    @CIRCachedProperty
    def rule_of_5_violation_count(self): return resolve(self.input, 'rule_of_5_violation_count', self.resolvers)

    @CIRCachedProperty
    def rotor_count(self): return resolve(self.input, 'rotor_count', self.resolvers)

    @CIRCachedProperty
    def effective_rotor_count(self): return resolve(self.input, 'effective_rotor_count', self.resolvers)

    @CIRCachedProperty
    def ring_count(self): return resolve(self.input, 'ring_count', self.resolvers)

    @CIRCachedProperty
    def ringsys_count(self): return resolve(self.input, 'ringsys_count', self.resolvers)



    @property
    def image_url(self):
        url = API_BASE+'/%s/image' % self.input
        if self.resolvers is not None:
            r = ",".join(self.resolvers)
            url += '?resolver=%s' % r
        return url

    @property
    def twirl_url(self):
        url = API_BASE+'/%s/twirl' % self.input
        if self.resolvers is not None:
            r = ",".join(self.resolvers)
            url += '?resolver=%s' % r
        return url

    @property
    def sdf_url(self):
        url = API_BASE+'/%s/sdf' % self.input
        if self.resolvers is not None:
            r = ",".join(self.resolvers)
            url += '?resolver=%s' % r
        return url


    def download(self, filename, format='sdf', overwrite=False):
        """ Download the resolved structure as a file """
        download(self.input, filename, format, overwrite)


# TODO add more formats from http://cactus.nci.nih.gov/blog/?p=68
file_formats = [
    'xyz',
    'pdb',
    ]

for name in file_formats:
    getter = lambda obj, nme: get_in_format(obj.input, nme)
    getter.__name__ = name
    getter.__doc__ = "A string containing the Molecule in the '" + name + "' file format."
    setattr(CIRMolecule, name, CIRCachedFileFormat(getter))

