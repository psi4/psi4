"""Module containing functions that, when passed a qcdb.WrappedDatabase instance
*dbinstance*, return an array of reaction names that are a subset of
dbinstance.hrxn.keys(). Since the full database is passed in, reactions
can be filtered by molecule characteristics, reaction names, existing
subsets, etc. The official name of the subset is specified by the function
docstring. Second line of docstring becomes tagl.

"""


def genset_MXuDD(dbinstance):
    """mxdd
    near-equilibrium systems also in mxdd

    """
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
    """HB-5min
    near-equilibrium systems also in hb

    """
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
    """MX-5min
    near-equilibrium systems also in mx

    """
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
    """DD-5min
    near-equilibrium systems also in dd

    """
    try:
        ssA = set(dbinstance.sset['dd'].keys())
    except KeyError:
        ssA = set()
    try:
        ssB = set(dbinstance.sset['5min'].keys())
    except KeyError:
        ssB = set()
    return ssA.intersection(ssB)


def genset_MXDDPPn5min(dbinstance):
    """MXDDPP-5min
    near-equilibrium systems also in mxddpp

    """
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
    """MXDDNP-5min
    near-equilibrium systems also in mxddnp

    """
    try:
        ssA = set(dbinstance.sset['mxddnp'].keys())
    except KeyError:
        ssA = set()
    try:
        ssB = set(dbinstance.sset['5min'].keys())
    except KeyError:
        ssB = set()
    return ssA.intersection(ssB)


def genset_allneutral(dbinstance):
    """neutral
    systems where all components are neutral

    """
    eligible = []
    for rxn, orxn in dbinstance.hrxn.iteritems():
        if all([True if rgt.charge == 0 else False for rgt in orxn.rxnm['default'].keys()]):
            eligible.append(rxn)
    return eligible


def genset_anyanion(dbinstance):
    """anion
    systems where any component is an anion

    """
    eligible = []
    for rxn, orxn in dbinstance.hrxn.iteritems():
        for rgt in orxn.rxnm['default'].keys():
            if rgt.charge < 0:
                eligible.append(rxn)
                break
    return eligible


def genset_anycation(dbinstance):
    """cation
    systems where any component is a cation

    """
    eligible = []
    for rxn, orxn in dbinstance.hrxn.iteritems():
        for rgt in orxn.rxnm['default'].keys():
            if rgt.charge > 0:
                eligible.append(rxn)
                break
    return eligible
