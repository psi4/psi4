"""Module containing functions that, when passed a qcdb.Database instance
*dbinstance*, return an array of reaction names that are a subset of
dbinstance.hrxn.keys(). Since the full database is passed in, reactions
can be filtered by molecule characteristics, reaction names, existing
subsets, etc. The official name of the subset is specified by the function
docstring.

"""


def genset_MXuDD(dbinstance):
    """mxdd"""
    try:
        ssA = set(dbinstance.sset['mx'].keys())
    except KeyError:
        ssA = set()
    try:
        ssB = set(dbinstance.sset['dd'].keys())
    except KeyError:
        ssB = set()
    return ssA.union(ssB)


def genset_HBn5min(dbinstance):
    """HB-5min"""
    try:
        ssA = set(dbinstance.sset['hb'].keys())
    except KeyError:
        ssA = set()
    try:
        ssB = set(dbinstance.sset['5min'].keys())
    except KeyError:
        ssB = set()
    return ssA.intersection(ssB)


def genset_MXn5min(dbinstance):
    """MX-5min"""
    try:
        ssA = set(dbinstance.sset['mx'].keys())
    except KeyError:
        ssA = set()
    try:
        ssB = set(dbinstance.sset['5min'].keys())
    except KeyError:
        ssB = set()
    return ssA.intersection(ssB)


def genset_DDn5min(dbinstance):
    """DD-5min"""
    try:
        ssA = set(dbinstance.sset['dd'].keys())
    except KeyError:
        ssA = set()
    try:
        ssB = set(dbinstance.sset['5min'].keys())
    except KeyError:
        ssB = set()
    return ssA.intersection(ssB)


#def genset_oldMXuDD(dbinstance):
#    """oldmxdd"""
#    try:
#        ssA = set(dbinstance.sset['oldmx'].keys())
#    except KeyError:
#        ssA = set()
#    try:
#        ssB = set(dbinstance.sset['oldhb'].keys())
#    except KeyError:
#        ssB = set()
#    return ssA.union(ssB)


def genset_MXDDPPn5min(dbinstance):
    """MXDDPP-5min"""
    try:
        ssA = set(dbinstance.sset['mxddpp'].keys())
    except KeyError:
        ssA = set()
    try:
        ssB = set(dbinstance.sset['5min'].keys())
    except KeyError:
        ssB = set()
    return ssA.intersection(ssB)


def genset_MXDDNPn5min(dbinstance):
    """MXDDNP-5min"""
    try:
        ssA = set(dbinstance.sset['mxddnp'].keys())
    except KeyError:
        ssA = set()
    try:
        ssB = set(dbinstance.sset['5min'].keys())
    except KeyError:
        ssB = set()
    return ssA.intersection(ssB)
