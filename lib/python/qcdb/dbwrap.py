import os
import sys
import math
try:
    import cPickle as pickle
except ImportError:
    import pickle
import itertools
from exceptions import *
from modelchems import Method, BasisSet, Error, methods, bases, errors
import psiutil
import textables


def initialize_errors(e=None, pe=None, pbe=None, extrema=True):
    """

    """
    error = OrderedDict()
    error['maxe'] = None if (e is None or not extrema) else e        # LD_XA
    error['mine'] = None if (e is None or not extrema) else e        # LD_XI
    error['me'] = None if e is None else 0.0                         # LD_MS
    error['mae'] = None if e is None else 0.0                        # LD_MA
    error['rmse'] = None if e is None else 0.0                       # LD_RA
    error['stde'] = None if e is None else 0.0
    error['maxpe'] = None if (pe is None or not extrema) else pe     # FD_XA
    error['minpe'] = None if (pe is None or not extrema) else pe     # FD_XI
    error['mpe'] = None if pe is None else 0.0                       # FD_MS
    error['mape'] = None if pe is None else 0.0                      # FD_MA
    error['rmspe'] = None if pe is None else 0.0                     # FD_RA
    error['stdpe'] = None if pe is None else 0.0
    error['maxpbe'] = None if (pbe is None or not extrema) else pbe  # BD_XA
    error['minpbe'] = None if (pbe is None or not extrema) else pbe  # BD_XI
    error['mpbe'] = None if pbe is None else 0.0                     # BD_MS
    error['mapbe'] = None if pbe is None else 0.0                    # BD_MA
    error['rmspbe'] = None if pbe is None else 0.0                   # BD_RA
    error['stdpbe'] = None if pbe is None else 0.0
    return error


def average_errors(*args):
    """Each item in *args* should be an error dictionary. Performs
    average-like operation over all items, which should be error
    dictionaries, in *args*. Defined for ME, MAE, STDE, and their
    relative-error variants. None returned for undefined statistics or
    when an item is missing.

    """
    Ndb = float(len(args))
    avgerror = initialize_errors()
    try:
        avgerror['maxe'] = max([x['maxe'] for x in args], key=lambda x: abs(x))
        avgerror['mine'] = min([x['mine'] for x in args], key=lambda x: abs(x))
        avgerror['me'] = sum([x['me'] for x in args]) / Ndb
        avgerror['mae'] = sum([x['mae'] for x in args]) / Ndb
        avgerror['rmse'] = 0.0  # TODO
        avgerror['stde'] = math.sqrt(sum([x['stde'] ** 2 for x in args]) / Ndb)

        avgerror['maxpe'] = max([x['maxpe'] for x in args], key=lambda x: abs(x))
        avgerror['minpe'] = min([x['minpe'] for x in args], key=lambda x: abs(x))
        avgerror['mpe'] = sum([x['mpe'] for x in args]) / Ndb
        avgerror['mape'] = sum([x['mape'] for x in args]) / Ndb
        avgerror['rmspe'] = 0.0  # TODO
        avgerror['stdpe'] = math.sqrt(sum([x['stdpe'] * x['stdpe'] for x in args]) / Ndb)

        avgerror['maxpbe'] = max([x['maxpbe'] for x in args], key=lambda x: abs(x))
        avgerror['minpbe'] = min([x['minpbe'] for x in args], key=lambda x: abs(x))
        avgerror['mpbe'] = sum([x['mpbe'] for x in args]) / Ndb
        avgerror['mapbe'] = sum([x['mapbe'] for x in args]) / Ndb
        avgerror['rmspbe'] = 0.0 # TODO
        avgerror['stdpbe'] = math.sqrt(sum([x['stdpbe'] * x['stdpbe'] for x in args]) / Ndb)
    except TypeError:
        pass
    return avgerror


def format_errors(err, mode=1):
    """From error dictionary *err*, returns a LaTeX-formatted string,
    after handling None entries.

    """
    if mode == 1:
        me = ' ----' if err['me'] is None else '%+.2f' % (err['me'])
        stde = '----' if err['stde'] is None else '%.2f' % (err['stde'])
        mae = '  ----' if err['mae'] is None else '%6.2f' % (err['mae'])
        mape = '  ----  ' if err['mape'] is None else '%6.1f\%%' % (100 * err['mape'])
        mapbe = '  ----  ' if err['mapbe'] is None else '%6.1f\%%' % (100 * err['mapbe'])
        text = """$\{%s; %s\}$ %s %s %s""" % \
            (me, stde, mae, mape, mapbe)
        return text

    if mode == 2:
        maxe = '----' if err['maxe'] is None else '%8.2f' % (err['maxe'])
        mine = '----' if err['mine'] is None else '%8.2f' % (err['mine'])
        me = '----' if err['me'] is None else '%+8.2f' % (err['me'])
        mae = '----' if err['mae'] is None else '%8.2f' % (err['mae'])
        rmse = '----' if err['rmse'] is None else '%8.2f' % (err['rmse'])
        stde = '----' if err['stde'] is None else '%8.2f' % (err['stde'])
        maxpe = '----' if err['maxpe'] is None else '%8.1f' % (100 * err['maxpe'])
        minpe = '----' if err['minpe'] is None else '%8.1f' % (100 * err['minpe'])
        mpe = '----' if err['mpe'] is None else '%+8.1f' % (100 * err['mpe'])
        mape = '----' if err['mape'] is None else '%8.1f' % (100 * err['mape'])
        rmspe = '----' if err['rmspe'] is None else '%8.1f' % (100 * err['rmspe'])
        stdpe = '----' if err['stdpe'] is None else '%8.1f' % (100 * err['stdpe'])
        maxpbe = '----' if err['maxpbe'] is None else '%8.1f' % (100 * err['maxpbe'])
        minpbe = '----' if err['minpbe'] is None else '%8.1f' % (100 * err['minpbe'])
        mpbe = '----' if err['mpbe'] is None else '%+8.1f' % (100 * err['mpbe'])
        mapbe = '----' if err['mapbe'] is None else '%8.1f' % (100 * err['mapbe'])
        rmspbe = '----' if err['rmspbe'] is None else '%8.1f' % (100 * err['rmspbe'])
        stdpbe = '----' if err['stdpbe'] is None else '%8.1f' % (100 * err['stdpbe'])
        text = """min: %s%s%s\nmax: %s%s%s\nm:   %s%s%s\nma:  %s%s%s\nrms: %s%s%s\nstd: %s%s%s""" % \
            (mine, minpe, minpbe, maxe, maxpe, maxpbe, me, mpe, mpbe, \
            mae, mape, mapbe, rmse, rmspe, rmspbe, stde, stdpe, stdpbe)
        return text

    if mode == 3:
        sdict = OrderedDict()
        sdict['maxe'] = '' if err['maxe'] is None else '%8.2f' % (err['maxe'])
        sdict['mine'] = '' if err['mine'] is None else '%8.2f' % (err['mine'])
        sdict['me'] = '' if err['me'] is None else '%+8.2f' % (err['me'])
        sdict['mae'] = '' if err['mae'] is None else '%8.2f' % (err['mae'])
        sdict['rmse'] = '' if err['rmse'] is None else '%8.2f' % (err['rmse'])
        sdict['stde'] = '' if err['stde'] is None else '%8.2f' % (err['stde'])
        sdict['maxpe'] = '' if err['maxpe'] is None else '%8.1f' % (100 * err['maxpe'])
        sdict['minpe'] = '' if err['minpe'] is None else '%8.1f' % (100 * err['minpe'])
        sdict['mpe'] = '' if err['mpe'] is None else '%+8.1f' % (100 * err['mpe'])
        sdict['mape'] = '' if err['mape'] is None else '%8.1f' % (100 * err['mape'])
        sdict['rmspe'] = '' if err['rmspe'] is None else '%8.1f' % (100 * err['rmspe'])
        sdict['stdpe'] = '' if err['stdpe'] is None else '%8.1f' % (100 * err['stdpe'])
        sdict['maxpbe'] = '' if err['maxpbe'] is None else '%8.1f' % (100 * err['maxpbe'])
        sdict['minpbe'] = '' if err['minpbe'] is None else '%8.1f' % (100 * err['minpbe'])
        sdict['mpbe'] = '' if err['mpbe'] is None else '%+8.1f' % (100 * err['mpbe'])
        sdict['mapbe'] = '' if err['mapbe'] is None else '%8.1f' % (100 * err['mapbe'])
        sdict['rmspbe'] = '' if err['rmspbe'] is None else '%8.1f' % (100 * err['rmspbe'])
        sdict['stdpbe'] = '' if err['stdpbe'] is None else '%8.1f' % (100 * err['stdpbe'])
        return sdict


def string_contrast(ss):
    """From an array of strings, *ss*, returns maximum common prefix
    string, maximum common suffix string, and array of middles.

    """
    s = [item + 'q' for item in ss if item is not None]
    short = min(s, key=len)
    for ib in range(len(short)):
        if not all([mc[ib] == short[ib] for mc in s]):
            preidx = ib
            break
    else:
        preidx = 0
    for ib in range(len(short)):
        ie = -1 * (ib + 1)
        if not all([mc[ie] == short[ie] for mc in s]):
            sufidx = ie + 1
            break
    else:
        sufidx = -1 * (len(short))

    miditer = iter([mc[preidx:sufidx] for mc in s])
    prefix = short[:preidx]
    suffix = short[sufidx:-1]
    middle = ['' if mc is None else next(miditer) for mc in ss]

    return prefix, suffix, middle


class ReactionDatum(object):
    """Piece of quantum chemical information that describes a qcdb.Reaction object.

    """
    def __init__(self, dbse, rxn, method, mode, basis, value, units='kcal/mol', comment=None):
        # geometry
        self.dbrxn = dbse + '-' + str(rxn)
        # qcdb.Method
        self.method = method
        # mode, e.g., unCP, CP, RLX, etc.
        self.mode = mode
        # qcdb.BasisSet
        self.basis = basis
        # numerical value for reaction
        self.value = float(value)
        # energy unit attached to value, defaults to kcal/mol
        self.units = units
        # addl comments
        self.comment = comment

    @classmethod
    def library_modelchem(cls, dbse, rxn, method, mode, basis, value, units='kcal/mol', comment=None):
        """Constructor when method and basis are strings corresponding to
        qcdb.Method and qcdb.BasisSet already defined in methods and bases.

        """
        # computational method
        if method.upper() in methods:
            tmp_method = methods[method.upper()]
        else:
            raise ValidationError("""Invalid ReactionDatum method %s.""" % (method))
        # computational basis set
        if basis.lower() in bases:
            tmp_basis = bases[basis.lower()]
        else:
            raise ValidationError("""Invalid ReactionDatum basis %s.""" % (basis))
        return cls(dbse, rxn, tmp_method, mode, tmp_basis, value, units='kcal/mol', comment=None)

    def __str__(self):
        text = ''
        text += """  ==> ReactionDatum <==\n\n"""
        text += """  Database reaction:    %s\n""" % (self.dbrxn)
        text += """  Method:               %s\n""" % (self.method.fullname)
        text += """  Mode:                 %s\n""" % (self.mode)
        text += """  Basis:                %s\n""" % (self.basis.fullname)
        text += """  Value:                %f [%s]\n""" % (self.value, self.units)
        text += """  Comment:              %s\n""" % (self.comment)
        text += """\n"""
        return text


class Reagent(object):
    """Chemical entity only slightly dresed up from qcdb.Molecule.

    """

    def __init__(self, name, mol, tagl=None, comment=None):
        # full name, e.g., 'S22-2-dimer' or 'NBC1-BzMe-8.0-monoA-CP' or 'HTBH-HCl-reagent'
        self.name = name
        # qcdb.Molecule
        try:
            self.NRE = mol.nuclear_repulsion_energy()
        except AttributeError:
            raise ValidationError("""Reagent must be instantiated with qcdb.Molecule object.""")
        else:
            self.mol = mol
        # description line
        self.tagl = tagl
        # addl comments
        self.comment = comment

    def __str__(self):
        text = ''
        text += """  ==> %s Reagent <==\n\n""" % (self.name)
        text += """  Tagline:              %s\n""" % (self.tagl)
        text += """  Comment:              %s\n""" % (self.comment)
        text += """  NRE:                  %f\n""" % (self.NRE)
        text += """  Molecule:             \n%s""" % (self.mol.format_molecule_for_psi4())
        text += """\n"""
        return text


class Reaction(object):
    """

    """

    def __init__(self, name, dbse, indx, tagl=None, latex=None, color='black', comment=None):
        # name, e.g., '2' or 'BzMe-8.0'
        self.name = name
        # database reaction name, e.g., 'S22-2' or 'NBC1-BzMe-8.0'
        self.dbrxn = dbse + '-' + str(name)
        # numerical index of reaction
        self.indx = indx
        # description line
        self.tagl = tagl
        # latex description
        self.latex = latex
        # addl comments
        self.comment = comment
        # reaction matrices, specifying reagent contributions per reaction
        self.rxnm = {}
        # qcdb.ReactionDatum objects of quantum chemical data pertaining to reaction
        self.data = {}
        # benchmark qcdb.ReactionDatum
        self.benchmark = None
        # color for plotting
        self.color = color

    def __str__(self):
        text = ''
        text += """  ==> %s Reaction <==\n\n""" % (self.name)
        text += """  Database reaction:    %s\n""" % (self.dbrxn)
        text += """  Index:                %s\n""" % (self.indx)
        text += """  LaTeX representation: %s\n""" % (self.latex)
        text += """  Tagline:              %s\n""" % (self.tagl)
        text += """  Comment:              %s\n""" % (self.comment)
        text += """  Benchmark:            %f\n""" % (self.data[self.benchmark].value)
        text += """  Color:                %s\n""" % (str(self.color))
        text += """  Reaction matrix:\n"""
        for mode, rxnm in self.rxnm.iteritems():
            text += """      %s\n""" % (mode)
            for rgt, coeff in rxnm.iteritems():
                text += """       %3d  %s\n""" % (coeff, rgt.name)
        text += """  Data:\n"""
        for label, datum in sorted(self.data.iteritems()):
            text += """      %8.2f  %s\n""" % (datum.value, label)
        text += """\n"""
        return text


class WrappedDatabase(object):
    """Wrapper class for raw Psi4 database modules that does some validation 
    of contents, creates member data and accessors for database structures, 
    defines error computation, and handles database subsets. Not to be used 
    directly-- see qcdb.Database for handling single or multiple 
    qdcb.WrappedDatabase objects and defining nice statistics, plotting, and
    table functionalities.

    >>> asdf = qcdb.WrappedDatabase('Nbc10')
    """

    def __init__(self, dbname, pythonpath=None):
        """Instantiate class with case insensitive name *dbname*. Module 
        search path can be prepended with *pythonpath*.

        """
        #: internal name of database
        #:
        #: >>> print asdf.dbse
        #: 'NBC1'
        self.dbse = None

        #: OrderedDict of reactions/members
        #:
        #: >>> print asdf.hrxn.keys()
        #: ['BzBz_S-3.2', 'BzBz_S-3.3', ...  'BzBz_PD36-2.8', 'BzBz_PD36-3.0']
        self.hrxn = None

        #: dict of reagents/geometries
        #:
        #: >>> print asdf.hrgt.keys()
        #: ['NBC1-BzBz_PD32-0.8-monoA-CP', 'NBC1-BzBz_PD34-0.6-dimer', ... 'NBC1-BzBz_PD34-1.7-dimer']
        self.hrgt = None

        #: dict of defined reaction subsets.
        #: Note that self.sset['default'] contains all the nonredundant information.
        #:
        #: >>> print asdf.sset.keys()
        #: ['meme', 'mxddpp', '5min', ... 'small']
        self.sset = None

        #   Removing hrxn, hrgt etc. do not reduce the size of the object.
        #   These attributes are stored for ease of access for adding qc info, etc.

        # load database
        if pythonpath is not None:
            sys.path.insert(1, pythonpath)
        else:
            sys.path.append(os.path.dirname(__file__) + '/../databases')
        try:
            database = psiutil.import_ignorecase(dbname)
        except ImportError:
            print('\nPython module for database %s failed to load\n\n' % (dbname))
            print('\nSearch path that was tried:\n')
            print(", ".join(map(str, sys.path)))
            raise ValidationError("Python module loading problem for database " + str(dbname))

        # gross validation of database
        for item in ['dbse', 'GEOS', 'HRXN', 'ACTV', 'RXNM']:
            try:
                getattr(database, item)
            except AttributeError:
                raise ValidationError("""Database %s severely deformed with %s missing.""" % (database.__name__, item))
        for item in ['TAGL', 'BIND']:
            try:
                getattr(database, item)
            except AttributeError:
                print """Warning: Database %s possibly deformed with %s missing.\n""" % (database.__name__, item)
        self.dbse = database.dbse

        # form array of database contents to process through
        pieces = []
        for item in dir(database):
            if item in ['qcdb', 'rxn', 'dbse', 'TAGL']:
                pass
            elif item.startswith('__'):
                pass
            else:
                pieces.append(item)

        # form qcdb.Reagent objects from all defined geometries, GEOS
        oHRGT = {}
#        for rgt, mol in database.GEOS.iteritems():
#            mol.update_geometry()
#            try:
#                tagl = database.TAGL[rgt]
#            except KeyError:
#                tagl = None
#                print """Warning: TAGL missing for reagent %s""" % (rgt)
#            oHRGT[rgt] = Reagent(name=rgt, mol=mol, tagl=tagl)
#        pieces.remove('GEOS')
        self.hrgt = oHRGT

        # form qcdb.Reaction objects from comprehensive reaction list, HRXN
        oHRXN = OrderedDict()
        for rxn in database.HRXN:
            try:
                tagl = database.TAGL[database.dbse + '-' + str(rxn)]
            except KeyError:
                tagl = None
                print """Warning: TAGL missing for reaction %s""" % (rxn)
            try:
                elst = database.DATA['SAPT ELST ENERGY'][database.dbse + '-' + str(rxn)]
                disp = database.DATA['SAPT DISP ENERGY'][database.dbse + '-' + str(rxn)]
                color = abs(elst) / (abs(elst) + abs(disp))
            except (KeyError, AttributeError):
                color = 'black'
                print """Warning: DATA['SAPT * ENERGY'] missing for reaction %s""" % (rxn)

            oHRXN[rxn] = Reaction(name=rxn,
                                  dbse=database.dbse,
                                  indx=database.HRXN.index(rxn) + 1,
                                  color=color,
                                  tagl=tagl)
        pieces.remove('HRXN')
        self.hrxn = oHRXN

        # list and align database stoichiometry modes, ACTV* and RXNM*
        oACTV = {}
        for modactv in [item for item in pieces if item.startswith('ACTV')]:
            modrxnm = modactv.replace('ACTV', 'RXNM')
            mode = 'default' if modactv == 'ACTV' else modactv.replace('ACTV_', '')
            try:
                getattr(database, modrxnm)
            except AttributeError:
                modrxnm = 'RXNM'
            oACTV[mode] = [modactv, modrxnm]
        for item in [tmp for tmp in pieces if tmp.startswith('ACTV') or tmp.startswith('RXNM')]:
            pieces.remove(item)

        # populate reaction matrices in qcdb.Reaction objects
#        for rxn in database.HRXN:
#            dbrxn = database.dbse + '-' + str(rxn)
#            for mode, actvrxnm in oACTV.iteritems():
#                tdict = OrderedDict()
#                for rgt in getattr(database, actvrxnm[0])[dbrxn]:
#                    tdict[oHRGT[rgt]] = getattr(database, actvrxnm[1])[dbrxn][rgt]
#                oHRXN[rxn].rxnm[mode] = tdict

        # list embedded quantum chem info per rxn, incl. BIND*
        arrsbind = [item for item in pieces if item.startswith('BIND_')]
        if len(arrsbind) == 0:
            if 'BIND' in pieces:
                arrsbind = ['BIND']
            else:
                arrsbind = []
                print """Warning: No BIND array with reference values."""
        else:
            for arrbind in arrsbind:
                if getattr(database, arrbind) is database.BIND:
                    break
            else:
                print """Warning: No BIND_* array assigned to be master BIND."""

        oBIND = {}
        for arrbind in arrsbind:
            ref = database.dbse + 'REF' if arrbind == 'BIND' else arrbind.replace('BIND_', '')
            methods[ref] = Method(name=ref)
            bases[ref] = BasisSet(name=ref)
            oBIND[ref] = [methods[ref], 'default', bases[ref], arrbind,
                (getattr(database, arrbind) is database.BIND)]
        for item in [tmp for tmp in pieces if tmp.startswith('BIND')]:
            pieces.remove(item)

        # populate data with reference values in qcdb.Reaction objects
        for rxn in database.HRXN:
            dbrxn = database.dbse + '-' + str(rxn)
            for ref, info in oBIND.iteritems():
                oHRXN[rxn].data[ref] = ReactionDatum(dbse=database.dbse,
                                                     rxn=rxn,
                                                     method=info[0],
                                                     mode=info[1],
                                                     basis=info[2],
                                                     value=getattr(database, info[3])[dbrxn])
                if info[4]:
                    oHRXN[rxn].benchmark = ref

        # Process subsets
        oSSET = {}
        fsHRXN = frozenset(database.HRXN)
        for sset in pieces:
            try:
                fssset = frozenset(getattr(database, sset))
            except TypeError:
                continue
            if fssset.issubset(fsHRXN):
                oSSET[sset] = getattr(database, sset)
        for item in oSSET.keys():
            pieces.remove(item)
        oSSET['HRXN'] = database.HRXN

        self.sset = OrderedDict()
        for item in oSSET.keys():
            if item == 'HRXN_SM':
                label = 'small'
            elif item == 'HRXN_LG':
                label = 'large'
            elif item == 'HRXN_EQ':
                label = 'equilibrium'
            elif item == 'HRXN':
                label = 'default'
            elif item.startswith('HRXN_'):
                label = item.replace('HRXN_', '').lower()
            else:
                label = item.lower()

            # subsets may have different ordering from HRXN
            self.sset[label] = OrderedDict()
            for rxn in oSSET[item]:
                self.sset[label][rxn] = oHRXN[rxn]

        print """WrappedDatabase %s: Unparsed attributes""" % (self.dbse), pieces

    def __str__(self):
        text = ''
        text += """  ==> %s WrappedDatabase <==\n\n""" % (self.dbse)
        text += """  Reagents:             %s\n""" % (self.hrgt.keys())
        text += """  Reactions:            %s\n""" % (self.hrxn.keys())
        text += """  Subsets:              %s\n""" % (self.sset.keys())
        text += """  Reference:            %s\n""" % (self.benchmark())
        text += """\n"""
        return text

    def add_ReactionDatum(self, dbse, rxn, method, mode, basis, value, units='kcal/mol', comment=None, overwrite=False):
        """Add a new quantum chemical value to *rxn* by creating a
        qcdb.ReactionDatum from same arguments as that class's
        object-less constructor. *rxn* may be actual Reaction.name
        or Reaction.indx.

        """
        if (self.dbse == dbse):
            if rxn in self.hrxn.keys():
                rxnname = rxn  # rxn is proper reaction name
            else:
                try:
                    if (rxn + 1 > 0) and (rxn == self.hrxn.items()[rxn - 1][1].indx):
                        rxnname = self.hrxn.items()[rxn - 1][1].name  # rxn is reaction index (maybe dangerous?)
                except (TypeError, IndexError):
                    raise ValidationError("""Inconsistent to add ReactionDatum for %s to database %s with reactions %s.""" %
                        (dbse + '-' + str(rxn), self.dbse, self.hrxn.keys()))
            label = '-'.join([method, mode, basis])
            if overwrite or (label not in self.hrxn[rxnname].data.keys()):
                self.hrxn[rxnname].data[label] = ReactionDatum.library_modelchem(dbse=dbse, rxn=rxnname,
                    method=method, mode=mode, basis=basis,
                    value=value, units=units, comment=comment)
            else:
                raise ValidationError("""ReactionDatum %s already present in Database.""" % (label))
        else:
            raise ValidationError("""Inconsistent to add ReactionDatum for %s to database %s.""" %
                (dbse + '-' + str(rxn), self.dbse))

    def add_Subset(self, name, func):
        """Define a new subset labeled *name* by providing a function
        *func* that filters *self.hrxn*.

        """
        label = name.lower()
        try:
            lsslist = [rxn for rxn in self.sset['default'].keys() if rxn in func(self)]
        except TypeError, e:
            raise ValidationError("""Function %s did not return list: %s.""" % (func.__name__, str(e)))
        if len(lsslist) == 0:
            print """WrappedDatabase %s: Subset %s NOT formed: empty""" % (self.dbse, label)
            return

        self.sset[label] = OrderedDict()
        for rxn in lsslist:
            self.sset[label][rxn] = self.hrxn[rxn]
        print """WrappedDatabase %s: Subset %s formed: %s""" % (self.dbse, label, self.sset[label].keys())

    def compute_errors(self, modelchem, benchmark='default', sset='default', failoninc=True, verbose=False):
        """For full database or subset *sset*, computes raw reaction
        errors between *modelchem* and *benchmark* model chemistries.
        Returns error if model chemistries are missing for any reaction in
        subset unless *failoninc* set to False, whereupon returns partial.
        Returns dictionary of reaction labels and error forms.

        """
        if isinstance(sset, basestring):
            # sset is normal subset name 'MX' corresponding to HRXN_MX or MX array in database module
            try:
                lsset = self.sset[sset.lower()]
            except KeyError, e:
                #raise ValidationError("""Subset named %s not available""" % (str(e)))
                lsset = OrderedDict()
        else:
            if callable(sset):
                # sset is function that will generate subset of HRXN from sset(self)
                lsslist = [rxn for rxn in self.sset['default'].keys() if rxn in sset(self)]
            else:
                # sset is array containing reactions
                lsslist = [rxn for rxn in self.sset['default'].keys() if rxn in sset]
            # assemble dict of qcdb.Reaction objects from array of reaction names
            lsset = OrderedDict()
            for rxn in lsslist:
                lsset[rxn] = self.hrxn[rxn]

        err = {}
        for rxn, oRxn in lsset.iteritems():
            lbench = oRxn.benchmark if benchmark == 'default' else benchmark
            try:
                mcLesser = oRxn.data[modelchem].value
                mcGreater = oRxn.data[lbench].value
            except KeyError, e:
                if failoninc:
                    raise ValidationError("""Reaction %s missing datum %s.""" % (str(rxn), str(e)))
                else:
                    continue

            err[rxn] = [mcLesser - mcGreater,
                        (mcLesser - mcGreater) / abs(mcGreater),
                        (mcLesser - mcGreater) / abs(mcGreater)]  # TODO define BER
            if verbose:
                print """p = %6.2f, pe = %6.1f%%, bpe = %6.1f%% reaction %s.""" % \
                    (err[rxn][0], 100 * err[rxn][1], 100 * err[rxn][2], str(rxn))
        return err

    def compute_statistics(self, modelchem, benchmark='default', sset='default', failoninc=True, verbose=False, returnindiv=False):
        """For full database or subset *sset*, computes many error
        statistics between single *modelchem* and *benchmark* model
        chemistries. Returns error if model chemistries are missing
        for any reaction in subset unless *failoninc* set to False,
        whereupon returns partial statistics. Returns dictionary of
        statistics labels and values.

        """
        err = self.compute_errors(modelchem, benchmark=benchmark, sset=sset, failoninc=failoninc, verbose=verbose)
        if len(err) == 0:
            error = initialize_errors()
            if verbose:
                print """Warning: nothing to compute."""
        else:
            Nrxn = float(len(err))
            error = OrderedDict()
            # linear (absolute) error
            linear = [val[0] for val in err.values()]
            error['maxe'] = max(linear, key=lambda x: abs(x))
            error['mine'] = min(linear, key=lambda x: abs(x))
            error['me'] = sum(linear) / Nrxn
            error['mae'] = sum(map(abs, linear)) / Nrxn
            error['rmse'] = math.sqrt(sum(map(lambda x: x ** 2, linear)) / Nrxn)
            error['stde'] = math.sqrt((sum(map(lambda x: x ** 2, linear)) - (sum(linear) ** 2) / Nrxn) / Nrxn)
            # fractional (relative) error
            relative = [val[1] for val in err.values()]
            error['maxpe'] = max(relative, key=lambda x: abs(x))
            error['minpe'] = min(relative, key=lambda x: abs(x))
            error['mpe'] = sum(relative) / Nrxn
            error['mape'] = sum(map(abs, relative)) / Nrxn
            error['rmspe'] = math.sqrt(sum(map(lambda x: x ** 2, relative)) / Nrxn)
            error['stdpe'] = math.sqrt((sum(map(lambda x: x ** 2, relative)) - (sum(relative) ** 2) / Nrxn) / Nrxn)
            # balanced (relative) error
            balanced = [val[2] for val in err.values()]
            error['maxpbe'] = max(balanced, key=lambda x: abs(x))
            error['minpbe'] = min(balanced, key=lambda x: abs(x))
            error['mpbe'] = sum(balanced) / Nrxn
            error['mapbe'] = sum(map(abs, balanced)) / Nrxn
            error['rmspbe'] = math.sqrt(sum(map(lambda x: x ** 2, balanced)) / Nrxn)
            error['stdpbe'] = math.sqrt((sum(map(lambda x: x ** 2, balanced)) - (sum(balanced) ** 2) / Nrxn) / Nrxn)
            if verbose:
                print """%d systems in %s for %s vs. %s, subset %s.\n%s""" % \
                    (len(err), self.dbse, modelchem, benchmark, sset, format_errors(error, mode=2))
        if returnindiv:
            return error, err
        else:
            return error

    def load_qcdata(self, modname, funcname, pythonpath=None, failoninc=True):
        """Loads qcdb.ReactionDatums from module *modname* function
        *funcname*. Module search path can be prepended with *pythonpath*.

        """
        if pythonpath is not None:
            sys.path.insert(1, pythonpath)
        else:
            sys.path.append(os.path.dirname(__file__) + '/../data')
        try:
            datamodule = __import__(modname)
        except ImportError:
            if not failoninc:
                print """%s data unavailable for database %s.\n""" % (modname, self.dbse)
                return
            else:
                print '\nPython module for database data %s failed to load\n\n' % (modname)
                print '\nSearch path that was tried:\n'
                print ", ".join(map(str, sys.path))
                raise ValidationError("Python module loading problem for database data " + str(modname))
        try:
            getattr(datamodule, funcname)(self)
        except AttributeError:
            if not failoninc:
                print """%s %s data unavailable for database %s.\n""" % (modname, funcname, self.dbse)
                return
            else:
                raise ValidationError("Python module missing function %s for loading data " % (str(funcname)))

        print """WrappedDatabase %s: %s %s results loaded""" % (self.dbse, modname, funcname)

    def load_qcdata_byproject(self, project, pythonpath=None):
        """Loads qcdb.ReactionDatums from standard location for *project*
        :module dbse_project and function load_project. Module search path
        can be prepended with *pythonpath*.

        """
        mod = self.dbse + '_' + project
        func = 'load_' + project
        self.load_qcdata(modname=mod, funcname=func, pythonpath=pythonpath)

    def load_qcdata_hdf5_trusted(self, project, path=None):
        """Loads qcdb.ReactionDatums from HDF5 file at path/dbse_project.h5 .
        If path not given, looks in qcdb/data. This file is written by
        reap-DB and so has been largely validated.

        """
        if path is None:
            path = os.path.dirname(__file__) + '/../data'
        hdf5file = os.path.abspath(path) + os.sep + self.dbse + '_' + project + '.h5'
        if not os.path.isfile(hdf5file):
            raise ValidationError("HDF5 file for loading database data from file %s does not exist" % (hdf5file))
        try:
            import pandas as pd
        except ImportError:
            raise ValidationError("Pandas data managment module must be available for import")

        try:
            self.hrxn.iterkeys().next() + 1
        except TypeError:
            intrxn = False
        else:
            intrxn = True

        with pd.get_store(hdf5file) as handle:
            for mc in handle['pdie'].keys():
                lmc = mc.split('-')  # TODO could be done better
                method = lmc[0]
                bsse = '_'.join(lmc[1:-1])
                basis = lmc[-1]

                df = handle['pdie'][mc]
                for dbrxn in df.index[df.notnull()].values:
                    [dbse, rxn] = dbrxn.split('-', 1)
                    if intrxn:
                        rxn = int(rxn)
                    self.hrxn[rxn].data[mc] = ReactionDatum.library_modelchem(dbse=dbse, rxn=rxn,
                        method=method, mode=bsse, basis=basis, value=df[dbrxn])

    @staticmethod
    def load_pickled(dbname, path=None):
        """

        """
        if path is None:
            path = os.path.dirname(__file__) + '/../data'
        picklefile = os.path.abspath(path) + os.sep + dbname + '.pickle'
        if not os.path.isfile(picklefile):
            raise ValidationError("Pickle file for loading database data from file %s does not exist" % (picklefile))
        with open(picklefile, 'rb') as handle:
            instance = pickle.load(handle)
        return instance

    def available_modelchems(self, union=True):
        """Returns all the labels of model chemistries that have been
        loaded. Either all modelchems that have data for any reaction if
        *union* is True or all modelchems that have data for all reactions
        if *union* is False.

        """
        mcs = [set(v.data) for k, v in self.hrxn.items()]
        if union:
            return sorted(set.union(*mcs))
        else:
            return sorted(set.intersection(*mcs))

    def benchmark(self):
        """Returns the model chemistry label for the database's benchmark."""
        return self.hrxn.itervalues().next().benchmark
        # TODO all rxns have same bench in db module so all have same here in obj
        #   but the way things stored in Reactions, this doesn't have to be so

    def load_subsets(self, modname='subsetgenerator', pythonpath=None):
        """Loads subsets from all functions in module *modname*.

        """
        if pythonpath is not None:
            sys.path.insert(1, pythonpath)
        else:
            sys.path.append(os.path.dirname(__file__))
        try:
            ssmod = __import__(modname)
        except ImportError:
            print '\nPython module for database data %s failed to load\n\n' % (modname)
            print '\nSearch path that was tried:\n'
            print ", ".join(map(str, sys.path))
            raise ValidationError("Python module loading problem for database subset generator " + str(modname))

        for func in dir(ssmod):
            if callable(getattr(ssmod, func)):
                self.add_Subset(getattr(ssmod, func).__doc__, getattr(ssmod, func))

        print """WrappedDatabase %s: Defined subsets loaded""" % (self.dbse)

    #def analyze_modelchems(self, modelchem, benchmark='default', failoninc=True, verbose=False):
    #    """Compute and print error statistics for each model chemistry in
    #    array *modelchem* versus *benchmark* for all available subsets and
    #    return dictionary of same.

    #    """
    #    pre, suf, mid = string_contrast(modelchem)
    #    errors = OrderedDict()
    #    for ss in self.sset.keys():
    #        errors[ss] = OrderedDict()
    #        for mc in modelchem:
    #            errors[ss][mc] = self.compute_statistics(mc, benchmark=benchmark, sset=ss, failoninc=failoninc, verbose=verbose)
    #    print """\n  ==> %s %s[]%s Errors <==""" % (self.dbse, pre, suf)
    #    print """%20s        %5s  %4s   %6s %6s    %6s""" % \
    #        ('', 'ME', 'STDE', 'MAE', 'MA%E', 'MA%BE')
    #    for ss in self.sset.keys():
    #        if any([any(errors[ss][mc].values()) for mc in modelchem]):
    #            print """   => %s <= """ % (ss)
    #            for mc in modelchem:
    #                print """%20s    %42s""" % (mid[modelchem.index(mc)], format_errors(errors[ss][mc]))
    #    return errors

    #def plot_modelchems(self, modelchem, benchmark='default', sset='default', failoninc=True, verbose=False, color='sapt', xlimit=4.0):
    #    """Computes individual errors and summary statistics for each model
    #    chemistry in array *modelchem* versus *benchmark* over subset *sset*.
    #    Thread *color* can be 'rgb' for old coloring, a color name or 'sapt'
    #    for spectrum coloring. Prepares thread diagram instructions and
    #    either executes them if matplotlib available (Canopy) or prints them.

    #    """
    #    pre, suf, mid = string_contrast(modelchem)
    #    # compute errors
    #    errors = OrderedDict()
    #    indiv = OrderedDict()
    #    for mc in modelchem:
    #        errors[mc], indiv[mc] = self.compute_statistics(mc, benchmark=benchmark,
    #            sset=sset, failoninc=failoninc, verbose=verbose, returnindiv=True)
    #    # repackage
    #    dbdat = []
    #    for rxn in self.sset[sset].keys():
    #        data = []
    #        for mc in modelchem:
    #            try:
    #                data.append(indiv[mc][rxn][0])
    #            except KeyError, e:
    #                if failoninc:
    #                    raise e
    #                else:
    #                    data.append(None)
    #        dbdat.append({'sys': str(rxn), 'color': self.hrxn[rxn].color, 'data': data})
    #    title = self.dbse + ' ' + pre + '[]' + suf
    #    mae = [errors[mc]['mae'] for mc in modelchem]
    #    mapbe = [100 * errors[mc]['mapbe'] for mc in modelchem]
    #    # generate matplotlib instructions and call or print
    #    try:
    #        import mpl
    #        import matplotlib.pyplot as plt
    #    except ImportError:
    #        # if not running from Canopy, print line to execute from Canopy
    #        print """mpl.thread(%s,\n    color='%s',\n    title='%s',\n    labels=%s,\n    mae=%s,\n    mape=%s\n    xlimit=%s)\n\n""" % \
    #            (dbdat, color, title, mid, mae, mapbe, str(xlimit))
    #    else:
    #        # if running from Canopy, call mpl directly
    #        mpl.thread(dbdat, color=color, title=title, labels=mid, mae=mae, mape=mapbe, xlimit=xlimit)

    #def plot_modelchems_mouseover(self, modelchem, benchmark='default', mbenchmark=None, sset='default', msset=None, failoninc=True, verbose=False, color='sapt', xlimit=4.0, saveas=None, mousetext=None, mouselink=None, mouseimag=None, mousetitle=None, force_relpath=False):
    #    """Computes individual errors and summary statistics for each model
    #    chemistry in array *modelchem* versus *benchmark* over subset *sset*.
    #    *mbenchmark* and *msset* are array options (same length as *modelchem*)
    #    that override *benchmark* and *sset*, respectively, for non-uniform
    #    specification. *saveas* conveys directory ('/') and/or filename
    #    for saving the resulting plot; file extension is not accessible.
    #    Thread *color* can be 'rgb' for old coloring, a color name or 'sapt'
    #    for spectrum coloring. Prepares thread diagram instructions and
    #    either executes them if matplotlib available (Canopy) or prints them.
    #    If any of *mousetext*, *mouselink*, or *mouseimag* is specified,
    #    htmlcode will be returned with an image map of slats to any of
    #    text, link, or image, respectively.

    #    """
    #    # distribute benchmark
    #    if mbenchmark is None:
    #        # benchmark is normal modelchem name
    #        lbenchmark = [benchmark] * len(modelchem)
    #    else:
    #        if isinstance(mbenchmark, basestring) or len(mbenchmark) != len(modelchem):
    #            raise ValidationError("""mbenchmark must be array of length distributable among modelchem""" % (str(mbenchmark)))
    #        else:
    #        # mbenchmark is array of benchmarks for each modelchem
    #            lbenchmark = mbenchmark
    #    # distribute sset
    #    if msset is None:
    #        # sset is normal subset name like 'MX'
    #        lsset = [sset] * len(modelchem)
    #    else:
    #        if isinstance(msset, basestring) or len(msset) != len(modelchem):
    #            raise ValidationError("""msset must be array of length distributable among modelchem""" % (str(msset)))
    #        else:
    #        # msset is array of subsets for each modelchem
    #            lsset = msset
    #    # compute errors
    #    errors = OrderedDict()
    #    indiv = OrderedDict()
    #    index = []
    #    for mc, bm, ss in zip(modelchem, lbenchmark, lsset):
    #        ix = '%s_%s_%s' % (mc, bm, ss)
    #        index.append(ix)
    #        errors[ix], indiv[ix] = self.compute_statistics(mc, benchmark=bm,
    #            sset=ss, failoninc=failoninc, verbose=verbose, returnindiv=True)
    #    # repackage
    #    dbdat = []
    #    for rxn in self.hrxn.keys():
    #        data = []
    #        for ix in index:
    #            if rxn in self.sset[lsset[index.index(ix)]].keys():
    #                try:
    #                    data.append(indiv[ix][rxn][0])
    #                except KeyError, e:
    #                    if failoninc:
    #                        raise e
    #                    else:
    #                        data.append(None)
    #            else:
    #                data.append(None)
    #        dbdat.append({'db': self.dbse, 'sys': str(rxn), 'color': self.hrxn[rxn].color, 'data': data})
    #    mae = [errors[ix]['mae'] for ix in index]
    #    mapbe = [100 * errors[ix]['mapbe'] for ix in index]
    #    # form unique filename
    #    ixpre, ixsuf, ixmid = string_contrast(index)
    #    title = self.dbse + ' ' + ixpre + '[]' + ixsuf
    #    # generate matplotlib instructions and call or print
    #    try:
    #        import mpl
    #        import matplotlib.pyplot as plt
    #    except ImportError:
    #        # if not running from Canopy, print line to execute from Canopy
    #        print """mpl.thread_mouseover_th(%s,\n    color='%s',\n    title='%s',\n    labels=%s,\n    mae=%s,\n    mape=%s\n    xlimit=%s\n    saveas=%s\n    mousetext=%s\n    mouselink=%s\n    mouseimag=%s\n    mousetitle=%s,\n    force_relpath=%s)\n\n""" % \
    #            (dbdat, color, title, ixmid, mae, mapbe, str(xlimit),
    #            repr(saveas), repr(mousetext), repr(mouselink), repr(mouseimag), repr(mousetitle), repr(force_relpath))
    #    else:
    #        # if running from Canopy, call mpl directly
    #        filedict, htmlcode = mpl.thread_mouseover(dbdat, color=color, title=title, labels=ixmid, mae=mae, mape=mapbe, xlimit=xlimit, saveas=saveas, mousetext=mousetext, mouselink=mouselink, mouseimag=mouseimag, mousetitle=mousetitle, force_relpath=force_relpath)
    #        return filedict, htmlcode

    #def plot_flat(self, modelchem, benchmark='default', sset='default', failoninc=True, verbose=False, color='sapt', xlimit=4.0, view=True):
    #    """Computes individual errors and summary statistics for single
    #    model chemistry *modelchem* versus *benchmark* over
    #    subset *sset*. Thread *color* can be 'rgb'
    #    for old coloring, a color name or 'sapt' for spectrum coloring.
    #    Prepares flat diagram instructions and either executes them if
    #    matplotlib available (Canopy) or prints them.

    #    """
    #    # compute errors
    #    errors = {}
    #    indiv = {}
    #    mc = modelchem
    #    errors[mc], indiv[mc] = self.compute_statistics(mc, benchmark=benchmark, sset=sset,
    #        failoninc=failoninc, verbose=verbose, returnindiv=True)
    #    # repackage
    #    dbdat = []
    #    for rxn in self.sset[sset].keys():
    #        data = []
    #        try:
    #            data.append(indiv[mc][rxn][0])
    #        except KeyError, e:
    #            if failoninc:
    #                raise e
    #            else:
    #                data.append(None)
    #        dbdat.append({'sys': str(rxn), 'color': self.hrxn[rxn].color, 'data': data})
    #    pre, suf, mid = string_contrast(mc)
    #    title = self.dbse + sset + ' ' + pre + '[]' + suf
    #    mae = errors[mc]['mae']
    #    mapbe = 100 * errors[mc]['mapbe']
    #    mapbe = None
    #    # generate matplotlib instructions and call or print
    #    try:
    #        import mpl
    #        import matplotlib.pyplot as plt
    #    except ImportError:
    #        # if not running from Canopy, print line to execute from Canopy
    #        print """mpl.flat(%s,\n    color='%s',\n    title='%s',\n    mae=%s,\n    mape=%s,\n    xlimit=%s,\n    view=%s)\n\n""" % \
    #            (dbdat, color, mc, mae, mapbe, xlimit, view)
    #    else:
    #        # if running from Canopy, call mpl directly
    #        mpl.flat(dbdat, color=color, title=mc, mae=mae, mape=mapbe, xlimit=xlimit, view=view)

    #def plot_bars(self, modelchem, benchmark='default', sset=['default', 'hb', 'mx', 'dd'], failoninc=True, verbose=False):
    #    """Prepares 'grey bars' diagram for each model chemistry in array
    #    *modelchem* versus *benchmark* over all four databases. A wide bar
    #    is plotted with three smaller bars, corresponding to the 'mae'
    #    summary statistic of the four subsets in *sset*. Prepares bars
    #    diagram instructions and either executes them if matplotlib
    #    available (Canopy) or prints them.

    #    """
    #    # compute errors
    #    errors = {}
    #    for mc in modelchem:
    #        if mc is not None:
    #            errors[mc] = {}
    #            for ss in sset:
    #                errors[mc][ss] = self.compute_statistics(mc, benchmark=benchmark, sset=ss,
    #                    failoninc=failoninc, verbose=verbose, returnindiv=False)
    #    # repackage
    #    pre, suf, mid = string_contrast(modelchem)
    #    #dbdat = [{'mc': mid[modelchem.index(mc)], 'data': [errors[mc][ss]['DB4']['mae'] for ss in sset]} for mc in modelchem]
    #    dbdat = []
    #    for mc in modelchem:
    #        if mc is None:
    #            dbdat.append(None)
    #        else:
    #            dbdat.append({'mc': mid[modelchem.index(mc)], 'data': [errors[mc][ss]['mae'] for ss in sset]})
    #    title = self.dbse + ' ' + pre + '[]' + suf
    #    # generate matplotlib instructions and call or print
    #    try:
    #        import mpl
    #        import matplotlib.pyplot as plt
    #    except ImportError:
    #        # if not running from Canopy, print line to execute from Canopy
    #        print """mpl.bar(%s,\n    title='%s')\n\n""" % (dbdat, title)
    #    else:
    #        # if running from Canopy, call mpl directly
    #        mpl.bar(dbdat, title=title)

    def plot_iowa(self, modelchem, benchmark='default', sset='default', failoninc=True, verbose=False, xlimit=2.0):
        """Computes individual errors for single *modelchem* versus
        *benchmark* over subset *sset*. Coloring green-to-purple with
        maximum intensity at *xlimit*. Prepares Iowa plot instructions and
        either executes them if matplotlib available (Canopy) or prints them.

        """
        title = self.dbse + ' ' + modelchem
        # compute errors
        errors, indiv = self.compute_statistics(modelchem, benchmark=benchmark,
                sset=sset, failoninc=failoninc, verbose=verbose, returnindiv=True)
        # repackage
        mcdat = []
        mclbl = []
        for rxn in self.sset[sset].keys():
            try:
                mcdat.append(indiv[rxn][0])
                mclbl.append(str(rxn))
            except KeyError, e:
                if failoninc:
                    raise e
        # generate matplotlib instructions and call or print
        try:
            import mpl
            import matplotlib.pyplot as plt
        except ImportError:
            # if not running from Canopy, print line to execute from Canopy
            print """mpl.iowa(%s,\n    %s,\n    title='%s',\n    xlimit=%s)\n\n""" % \
                (mcdat, mclbl, title, str(xlimit))
        else:
            # if running from Canopy, call mpl directly
            mpl.iowa(mcdat, mclbl, title=title, xlimit=xlimit)
            #print """mpl.iowa(%s,\n    %s,\n    title='%s',\n    xlimit=%s)\n\n""" % \
            #    (mcdat, mclbl, title, str(xlimit))

    def table_generic(self, mtd, bas, columnplan, rowplan=['bas', 'mtd'],
        opt=['CP'], err=['mae'], sset=['default'],
        benchmark='default', failoninc=True,
        landscape=False, standalone=True, subjoin=True,
        plotpath='', theme='', filename=None):
        """Prepares dictionary of errors for all combinations of *mtd*, *opt*,
        *bas* with respect to model chemistry *benchmark*, mindful of *failoninc*.
        Once error dictionary is ready, it and all other arguments are passed
        along to textables.table_generic.

        """
        # gather list of model chemistries for table
        mcs = ['-'.join(prod) for prod in itertools.product(mtd, opt, bas)]

        if plotpath == 'autogen':
            plotpath = os.environ['HOME'] + os.sep + 'mplflat_'
            for mc in mcs:
                self.plot_flat(mc)
            # TODO isn't going to work if sset in rowplan

        # compute errors
        serrors = {}
        for mc in mcs:
            serrors[mc] = {}
            for ss in self.sset.keys():
                errblock = self.compute_statistics(mc, benchmark=benchmark, sset=ss,
                    failoninc=failoninc, verbose=False, returnindiv=False)
                serrors[mc][ss] = {}
                serrors[mc][ss][self.dbse] = format_errors(errblock, mode=3)

        textables.table_generic(dbse=[self.dbse], serrors=serrors,
            mtd=mtd, bas=bas, columnplan=columnplan, rowplan=rowplan,
            opt=opt, err=err, sset=sset,
            landscape=landscape, standalone=standalone, subjoin=subjoin,
            plotpath=plotpath, theme=theme, filename=filename)

    def table_simple1(self, mtd, bas, opt=['CP'], err=['mae'], benchmark='default', failoninc=True, plotpath='analysis/flats/mplflat_', theme='smmerge'):
        """Specialization of table_generic into table with minimal statistics
        (three S22 and three overall) plus embedded slat diagram as suitable
        for main paper. A single table is formed in sections by *bas* with
        lines *mtd* within each section.

        """
        rowplan = ['bas', 'mtd']
        columnplan = [
            ['l', r"""Method \& Basis Set""", '', textables.label, {}],
            ['d', r'S22', 'HB', textables.val, {'sset': 'hb'}],
            ['d', r'S22', 'MX', textables.val, {'sset': 'mx'}],
            ['d', r'S22', 'DD', textables.val, {'sset': 'dd'}],
            ['d', r'S22', 'TT', textables.val, {'sset': 'default'}],
            ]

        self.table_generic(mtd=mtd, bas=bas, columnplan=columnplan, rowplan=rowplan,
            opt=opt, err=err,
            benchmark=benchmark, failoninc=failoninc,
            landscape=False, standalone=True, subjoin=True,
            plotpath=plotpath, theme=theme, filename=None)

    def table_simple2(self, mtd, bas, opt=['CP'], err=['mae'], benchmark='default', failoninc=True, plotpath='analysis/flats/mplflat_', theme='smmerge'):
        """Specialization of table_generic into table with minimal statistics
        (three S22 and three overall) plus embedded slat diagram as suitable
        for main paper. A single table is formed in sections by *bas* with
        lines *mtd* within each section.

        """
        rowplan = ['bas', 'mtd']
        columnplan = [
            ['l', r"""Method \& Basis Set""", '', textables.label, {}],
            ['d', r'MAE', 'HB', textables.val, {'sset': 'hb'}],
            ['d', r'MAE', 'MX', textables.val, {'sset': 'mx'}],
            ['d', r'MAE', 'DD', textables.val, {'sset': 'dd'}],
            ['d', r'MAE', 'TT', textables.val, {'sset': 'default'}],
            ['d', r'MA\%E', 'HB', textables.val, {'sset': 'hb', 'err': 'mape'}],
            ['d', r'MA\%E', 'MX', textables.val, {'sset': 'mx', 'err': 'mape'}],
            ['d', r'MA\%E', 'DD', textables.val, {'sset': 'dd', 'err': 'mape'}],
            ['d', r'MA\%E', 'TT', textables.val, {'sset': 'default', 'err': 'mape'}],
            ['d', r'maxE', 'TT ', textables.val, {'sset': 'default', 'err': 'maxe'}],
            ['d', r'min\%E', ' TT', textables.val, {'sset': 'default', 'err': 'minpe'}],
            ['d', r'rmsE', 'TT ', textables.val, {'sset': 'default', 'err': 'rmse'}],
            ['d', r'devE', ' TT', textables.val, {'sset': 'default', 'err': 'stde'}],
            ]

        self.table_generic(mtd=mtd, bas=bas, columnplan=columnplan, rowplan=rowplan,
            opt=opt, err=err,
            benchmark=benchmark, failoninc=failoninc,
            landscape=False, standalone=True, subjoin=True,
            plotpath=plotpath, theme=theme, filename=None)

    def table_simple3(self, mtd, bas, opt=['CP'], err=['mae'], benchmark='default', failoninc=True, plotpath='analysis/flats/mplflat_', theme='smmerge'):
        """Specialization of table_generic into table with minimal statistics
        (three S22 and three overall) plus embedded slat diagram as suitable
        for main paper. A single table is formed in sections by *bas* with
        lines *mtd* within each section.

        """
        rowplan = ['err', 'bas', 'mtd']
        columnplan = [
            ['l', r"""Method \& Basis Set""", '', textables.label, {}],
            ['d', r'MAE', 'HB', textables.val, {'sset': 'hb'}],
            ['d', r'MAE', 'MX', textables.val, {'sset': 'mx'}],
            ['d', r'MAE', 'DD', textables.val, {'sset': 'dd'}],
            ['d', r'MAE', 'TT', textables.val, {'sset': 'default'}],
            ]

        self.table_generic(mtd=mtd, bas=bas, columnplan=columnplan, rowplan=rowplan,
            opt=opt, err=err,
            benchmark=benchmark, failoninc=failoninc,
            landscape=False, standalone=True, subjoin=True,
            plotpath=plotpath, theme=theme, filename=None)

    def table_simple4(self, mtd, bas, opt=['CP'], err=['mae'], benchmark='default', failoninc=True, plotpath='analysis/flats/mplflat_', theme='smmerge'):
        """Specialization of table_generic into table with minimal statistics
        (three S22 and three overall) plus embedded slat diagram as suitable
        for main paper. A single table is formed in sections by *bas* with
        lines *mtd* within each section.

        """
        plotpath = 'autogen'  # TODO handle better
        rowplan = ['bas', 'mtd']
        columnplan = [
            ['l', r"""Method \& Basis Set""", '', textables.label, {}],
            ['d', r'S22', 'HB', textables.val, {'sset': 'hb'}],
            ['d', r'S22', 'MX', textables.val, {'sset': 'mx'}],
            ['d', r'S22', 'DD', textables.val, {'sset': 'dd'}],
            ['d', r'S22', 'TT', textables.val, {'sset': 'default'}],
            #['l', r"""Error Distribution\footnotemark[1]""", r"""\includegraphics[width=6.67cm,height=3.5mm]{%s%s.pdf}""" % (plotpath, 'blank'), textables.graphics, {}],
            ['l', r"""Error Distribution\footnotemark[1]""", r"""""", textables.graphics, {}],
            ]

        self.table_generic(mtd=mtd, bas=bas, columnplan=columnplan, rowplan=rowplan,
            opt=opt, err=err,
            benchmark=benchmark, failoninc=failoninc,
            landscape=False, standalone=True, subjoin=True,
            plotpath=plotpath, theme=theme, filename=None)


class Database(object):
    """Collection for handling single or multiple qcdb.WrappedDatabase objects.
    Particularly, unifying modelchem and subset names that when inconsistent
    across component databases. Also, defining statistics across databases.

    >>> asdf = qcdb.Database(['s22', 'Nbc10', 'hbc6', 'HSG'], 'DB4')
    >>> qwer = qcdb.Database(['s22'])
    """

    def __init__(self, dbnamelist, dbse=None, pythonpath=None, loadfrompickle=False, path=None):
        #: internal name of database collection
        #:
        #: >>> print asdf.dbse
        #: 'DB4'
        self.dbse = None

        #: ordered component Database objects
        #:
        #: >>> print asdf.dbdict
        #: XXXX
        self.dbdict = OrderedDict()

        #: subset assembly pattern
        #:
        #: >>> print asdf.sset.keys()
        #: XXXX
        self.sset = OrderedDict()

        #: assembly pattern for transspecies modelchems
        #:
        #: >>> print asdf.mcs.keys()
        #: XXXX
        self.mcs = {}

        self.benchmark = None

        # load databases
        for db in dbnamelist:
            if loadfrompickle:
                tmp = WrappedDatabase.load_pickled(db, path=path)
            else:
                tmp = WrappedDatabase(db, pythonpath=pythonpath)
            self.dbdict[tmp.dbse] = tmp

        # slurp up the obvious overlaps
        consolidated_bench = [odb.benchmark() for odb in self.dbdict.values()]
        if len(set(consolidated_bench)) == 1:
            self.benchmark = consolidated_bench[0]
        else:
            self.benchmark = ''.join(consolidated_bench)
        self.mcs[self.benchmark] = consolidated_bench

        #methods[ref] = Method(name=ref)
        #bases[ref] = BasisSet(name=ref)

        self.mcs['default'] = consolidated_bench
        #self.mcs['default'] = [odb.benchmark() for odb in self.dbdict.values()]
        self._intersect_subsets()
        self._intersect_modelchems()

        # collection name
        self.dbse = ''.join(self.dbdict.keys()) if dbse is None else dbse

        print """Database %s: %s""" % (self.dbse, ', '.join(self.dbdict.keys()))

    def __str__(self):
        text = ''
        text += """  ===> %s Database <===\n\n""" % (self.dbse)
        #text += """  Reagents:             %s\n""" % (self.hrgt.keys())
        #text += """  Reactions:            %s\n""" % (self.hrxn.keys())
        text += """  Subsets:              %s\n""" % (self.sset.keys())
        #text += """  Reference:            %s\n""" % ('default: ' + ' + '.join(self.mcs['default']))
        text += """  Reference:            %s\n""" % (self.benchmark + ': ' + ' + '.join(self.mcs[self.benchmark]))
        text += """  Model Chemistries:    %s\n""" % (', '.join(sorted(self.mcs.keys())))
        text += """\n"""
        for db in self.dbdict.keys():
            text += self.dbdict[db].__str__()
        return text

#    def benchmark(self):
#        """Returns the model chemistry label for the database's benchmark."""
#        return self.benchmark  #TODO not sure if right way to go about this self.mcs['default']

    def fancy_mcs(self):
        """

        """
        fmcs = {}
        for mc in self.mcs.keys():
            #print '%30s' % (mc),
            try:
                mtd, mod, bas = mc.split('-')
            except ValueError:
                fmcs[mc] = mc
            else:
                fmcs[mc] = """%20s / %-20s %s""" % (methods[mtd].fullname, bases[bas].fullname, mod)
            #print fmcs[mc]
        return fmcs

    def load_qcdata_byproject(self, project, pythonpath=None):
        """For each component database, loads qcdb.ReactionDatums from
        standard location for *project* :module dbse_project and function
        load_project. Module search path can be prepended with *pythonpath*.

        """
        for db, odb in self.dbdict.items():
            odb.load_qcdata_byproject(project, pythonpath=pythonpath)
        self._intersect_modelchems()

    def load_qcdata_hdf5_trusted(self, project, path=None):
        """For each component database, loads qcdb.ReactionDatums from
        HDF5 file at path/dbse_project.h5 . If path not given, looks in
        qcdb/data. This file is written by reap-DB and so has been largely
        validated.

        """
        for db, odb in self.dbdict.items():
            odb.load_qcdata_hdf5_trusted(project, path=path)
        self._intersect_modelchems()

    def load_subsets(self, modname='subsetgenerator', pythonpath=None):
        """For each component database, loads subsets from all functions
        in module *modname*. Default *modname* usues standard generators.

        """
        for db, odb in self.dbdict.items():
            odb.load_subsets(modname=modname, pythonpath=pythonpath)
        self._intersect_subsets()

    def _intersect_subsets(self):
        """Examine component database subsets and collect common names as
        Database subset.

        """
        sss = [set(odb.sset.keys()) for db, odb in self.dbdict.items()]
        new = sorted(set.intersection(*sss))
        for ss in new:
            self.sset[ss] = [ss] * len(self.dbdict.keys())

    def _intersect_modelchems(self):
        """Examine component database qcdata and collect common names as
        Database modelchem.

        """
        mcs = [set(odb.available_modelchems()) for db, odb in self.dbdict.items()]
        new = sorted(set.intersection(*mcs))
        for mc in new:
            self.mcs[mc] = [mc] * len(self.dbdict.keys())

    def compute_statistics(self, modelchem, benchmark='default', sset='default', failoninc=True, verbose=False, returnindiv=False):
        """Computes summary statistics and, if *returnindiv* True,
        individual errors for single model chemistry *modelchem* versus
        *benchmark* over subset *sset* over all component databases.
        Particularly, imposes cross-database definitions for sset and
        modelchem.
        #Returns error if model chemistries are missing
        #for any reaction in subset unless *failoninc* set to False,
        #whereupon returns partial statistics. Returns dictionary of
        #statistics labels and values.

        """
        errors = OrderedDict()
        indiv = OrderedDict()
        actvdb = []
        for db, odb in self.dbdict.items():
            dbix = self.dbdict.keys().index(db)
            if self.sset[sset][dbix] is None:
                errors[db], indiv[db] = (None, None)
            else:
                errors[db], indiv[db] = odb.compute_statistics(self.mcs[modelchem][dbix],
                    sset=self.sset[sset][dbix], benchmark=self.mcs[benchmark][dbix],
                    failoninc=failoninc, verbose=verbose, returnindiv=True)
                actvdb.append(errors[db])
        errors[self.dbse] = average_errors(*actvdb)

        if returnindiv:
            return errors, indiv
        else:
            return errors

    def analyze_modelchems(self, modelchem, benchmark='default', failoninc=True, verbose=False):
        """For each component database, compute and print nicely formatted
        summary error statistics for each model chemistry in array
        *modelchem* versus *benchmark* for all available subsets.

        """
        # compute errors
        errors = {}
        for mc in modelchem:
            errors[mc] = {}
            for ss in self.sset.keys():
                errors[mc][ss] = self.compute_statistics(mc, benchmark=benchmark, sset=ss,
                    failoninc=failoninc, verbose=verbose, returnindiv=False)
        # present errors
        pre, suf, mid = string_contrast(modelchem)
        text = """\n  ==> %s %s[]%s Errors <==\n""" % (self.dbse, pre, suf)
        text += """%20s    %44s""" % ('', '==> ' + self.dbse + ' <==')
        for db, odb in self.dbdict.items():
            text += """%44s""" % ('=> ' + odb.dbse + ' <=')
        text += '\n'
        text += """%20s         %5s   %4s   %6s %6s    %6s\n""" % \
            ('', 'ME', 'STDE', 'MAE', 'MA%E', 'MA%BE')
        for ss in self.sset.keys():
            text += """   => %s <=\n""" % (ss)
            for mc in modelchem:
                perr = errors[mc][ss]
                text += """%20s    %44s""" % (mid[modelchem.index(mc)],
                    format_errors(perr[self.dbse]))
                for db in self.dbdict.keys():
                    text += """%44s""" % ('' if perr[db] is None else format_errors(perr[db]))
                text += '\n'
        print text

    def plot_bars(self, modelchem, benchmark='default', sset=['default', 'hb', 'mx', 'dd'],
        failoninc=True, verbose=False,
        saveas=None, relpath=False, graphicsformat=['pdf']):
        """Prepares 'grey bars' diagram for each model chemistry in array
        *modelchem* versus *benchmark* over all component databases. A wide bar
        is plotted with three smaller bars, corresponding to the 'mae'
        summary statistic of the four subsets in *sset*.

        *saveas* conveys directory ('/') and/or filename for saving the
        resulting plot. File extension is not accessible, but *graphicsformat*
        array requests among 'png', 'pdf', and 'eps' formats. *relpath*
        forces paths to saved files to be relative to current directory,
        rather than absolute paths for returned code and file dictionary.

        Prepares bars diagram instructions and either executes them if
        matplotlib available (Canopy or Anaconda) or prints them. Returns a
        dictionary of all saved plot filenames.

        >>>> asdf.plot_bars(['MP2-CP-adz', 'MP2-CP-adtz'], sset=['tt-5min', 'hb-5min', 'mx-5min', 'dd-5min'])
        """
        # compute errors
        errors = {}
        for mc in modelchem:
            if mc is not None:
                errors[mc] = {}
                for ss in sset:
                    errors[mc][ss] = self.compute_statistics(mc, benchmark=benchmark, sset=ss,
                        failoninc=failoninc, verbose=verbose, returnindiv=False)
        # repackage
        pre, suf, mid = string_contrast(modelchem)
        dbdat = []
        for mc in modelchem:
            if mc is None:
                dbdat.append(None)
            else:
                dbdat.append({'mc': mid[modelchem.index(mc)],
                              'data': [errors[mc][ss][self.dbse]['mae'] for ss in sset]})
        title = self.dbse + ' ' + pre + '[]' + suf + ' ' + ','.join(sset)
        # generate matplotlib instructions and call or print
        try:
            import mpl
            import matplotlib.pyplot as plt
        except ImportError:
            # if not running from Canopy, print line to execute from Canopy
            print """filedict = mpl.bars(%s,\n    title='%s'\n    saveas=%s\n    relpath=%s\n    graphicsformat=%s)\n\n""" % \
                (dbdat, title, repr(saveas), repr(relpath), repr(graphicsformat))
        else:
            # if running from Canopy, call mpl directly
            filedict = mpl.bars(dbdat, title=title,
                saveas=saveas, relpath=relpath, graphicsformat=graphicsformat)
            return filedict

    def plot_flat(self, modelchem, benchmark='default', sset='default',
        failoninc=True, verbose=False, color='sapt', xlimit=4.0, view=True,
        saveas=None, relpath=False, graphicsformat=['pdf']):
        """Computes individual errors and summary statistics for single
        model chemistry *modelchem* versus *benchmark* over
        subset *sset* over all component databases. Thread *color* can be
        'rgb' for old coloring, a color name or 'sapt' for spectrum coloring.

        *saveas* conveys directory ('/') and/or filename for saving the
        resulting plot. File extension is not accessible, but *graphicsformat*
        array requests among 'png', 'pdf', and 'eps' formats. *relpath*
        forces paths to saved files to be relative to current directory,
        rather than absolute paths for returned code and file dictionary.

        Prepares flat diagram instructions and either executes them if
        matplotlib available (Canopy or Anaconda) or prints them. Returns a
        dictionary of all saved plot filenames.

        asdf.plot_flat('CCSD-CP-atqzadz', failoninc=False)
        """
        # compute errors
        mc = modelchem
        errors, indiv = self.compute_statistics(mc, benchmark=benchmark, sset=sset,
            failoninc=failoninc, verbose=verbose, returnindiv=True)
        # repackage
        dbdat = []
        for db, odb in self.dbdict.items():
            if indiv[db] is not None:
                for rxn in indiv[db].keys():
                    dbdat.append({'db': db,
                                  'sys': str(rxn),
                                  'color': odb.hrxn[rxn].color,
                                  'data': [indiv[db][rxn][0]]})
        pre, suf, mid = string_contrast(mc)
        title = self.dbse + '-' + sset + ' ' + pre + '[]' + suf
        mae = errors[self.dbse]['mae']
        mape = 100 * errors[self.dbse]['mape']
        mapbe = None
        # generate matplotlib instructions and call or print
        try:
            import mpl
            import matplotlib.pyplot as plt
        except ImportError:
            # if not running from Canopy, print line to execute from Canopy
            print """filedict = mpl.flat(%s,\n    color='%s',\n    title='%s',\n    mae=%s,\n    mape=%s,\n    xlimit=%s,\n    view=%s\n    saveas=%s\n    relpath=%s\n    graphicsformat=%s)\n\n""" % \
                (dbdat, color, mc, mae, mape, xlimit, view, repr(saveas), repr(relpath), repr(graphicsformat))
        else:
            # if running from Canopy, call mpl directly
            filedict = mpl.flat(dbdat, color=color, title=mc, mae=mae, mape=mape,
                xlimit=xlimit, view=view,
                saveas=saveas, relpath=relpath, graphicsformat=graphicsformat)
            return filedict

    def plot_all_flats(self, modelchem=None, sset='default', xlimit=4.0,
        saveas=None, relpath=False, graphicsformat=['pdf']):
        """Generate pieces for inclusion into tables. Supply list of
        modelchemistries to plot from *modelchem*, otherwise defaults to
        all those available. Can modify subset *sset* and plotting
        range *xlimit*.

        >>>> asdf.plot_all_flats(sset='tt-5min', xlimit=4.0)
        """
        mcs = self.mcs.keys() if modelchem is None else modelchem
        filedict = OrderedDict()
        for mc in sorted(mcs):
            minifiledict = self.plot_flat(mc, sset=sset, xlimit=xlimit, view=False,
                saveas=saveas, relpath=relpath, graphicsformat=graphicsformat)
            filedict[mc] = minifiledict['pdf']
        return filedict

    def plot_disthist(self, modelchem, benchmark='default', sset='default',
        failoninc=True, verbose=False, xtitle='',
        saveas=None, relpath=False, graphicsformat=['pdf']):
        """Computes individual errors and summary statistics for single
        model chemistry *modelchem* versus *benchmark* over
        subset *sset* over all component databases. Computes histogram
        of errors and gaussian distribution.

        *saveas* conveys directory ('/') and/or filename for saving the
        resulting plot. File extension is not accessible, but *graphicsformat*
        array requests among 'png', 'pdf', and 'eps' formats. *relpath*
        forces paths to saved files to be relative to current directory,
        rather than absolute paths for returned code and file dictionary.

        Prepares disthist diagram instructions and either executes them if
        matplotlib available (Canopy or Anaconda) or prints them. Returns a
        dictionary of all saved plot filenames.

        >>>>
        """
        # compute errors
        mc = modelchem
        errors, indiv = self.compute_statistics(mc, benchmark=benchmark, sset=sset,
            failoninc=failoninc, verbose=verbose, returnindiv=True)
        # repackage
        dbdat = []
        for db in self.dbdict.keys():
            if indiv[db] is not None:
                for rxn in indiv[db].keys():
                    dbdat.append(indiv[db][rxn][0])
        title = """%s vs %s for %s subset %s""" % (mc, benchmark, self.dbse, sset)
        me = errors[self.dbse]['me']
        stde = errors[self.dbse]['stde']
        # generate matplotlib instructions and call or print
        try:
            import mpl
            import matplotlib.pyplot as plt
        except ImportError:
            # if not running from Canopy, print line to execute from Canopy
            print """filedict = mpl.disthist(%s,\n    title='%s',\n    xtitle='%s'\n    me=%s,\n    stde=%s,\n    saveas=%s,\n    relpath=%s\n    graphicsformat=%s)\n\n""" % \
                (dbdat, title, xtitle, me, stde, repr(saveas), repr(relpath), repr(graphicsformat))
        else:
            # if running from Canopy, call mpl directly
            filedict = mpl.disthist(dbdat, title=title, xtitle=xtitle, me=me, stde=stde,
                saveas=saveas, relpath=relpath, graphicsformat=graphicsformat)
            return filedict

    def plot_modelchems(self, modelchem, benchmark='default', mbenchmark=None,
        sset='default', msset=None, failoninc=True, verbose=False, color='sapt',
        xlimit=4.0, saveas=None, mousetext=None, mouselink=None, mouseimag=None,
        mousetitle=None, relpath=False, graphicsformat=['pdf']):
        """Computes individual errors and summary statistics over all component
        databases for each model chemistry in array *modelchem* versus *benchmark*
        over subset *sset*. *mbenchmark* and *msset* are array options (same
        length as *modelchem*) that override *benchmark* and *sset*, respectively,
        for non-uniform specification. Thread *color* can be 'rgb' for old
        coloring, a color name or 'sapt' for spectrum coloring.

        *saveas* conveys directory ('/') and/or filename for saving the
        resulting plot. File extension is not accessible, but *graphicsformat*
        array requests among 'png', 'pdf', and 'eps' formats. *relpath*
        forces paths to saved files to be relative to current directory,
        rather than absolute paths for returned code and file dictionary.

        Prepares thread diagram instructions and either executes them if
        matplotlib available (Canopy or Anaconda) or prints them. Returns a
        dictionary of all saved plot filenames. If any of *mousetext*, *mouselink*,
        or *mouseimag* is specified, htmlcode will be returned with an image map of
        slats to any of text, link, or image, respectively.

        """
        # distribute benchmark
        if mbenchmark is None:
            lbenchmark = [benchmark] * len(modelchem)  # normal bm modelchem name
        else:
            if isinstance(mbenchmark, basestring) or len(mbenchmark) != len(modelchem):
                raise ValidationError("""mbenchmark must be array of length distributable among modelchem""" % (str(mbenchmark)))
            else:
                lbenchmark = mbenchmark  # array of bm for each modelchem
        # distribute sset
        if msset is None:
            lsset = [sset] * len(modelchem)  # normal ss name like 'MX'
        else:
            if isinstance(msset, basestring) or len(msset) != len(modelchem):
                raise ValidationError("""msset must be array of length distributable among modelchem""" % (str(msset)))
            else:
                lsset = msset  # array of ss for each modelchem
        # compute errors
        index = []
        errors = {}
        indiv = {}
        for mc, bm, ss in zip(modelchem, lbenchmark, lsset):
            ix = '%s_%s_%s' % (ss, mc, bm)
            index.append(ix)
            errors[ix], indiv[ix] = self.compute_statistics(mc, benchmark=bm, sset=ss,
                failoninc=failoninc, verbose=verbose, returnindiv=True)
        # repackage
        dbdat = []
        for db, odb in self.dbdict.items():
            for rxn in odb.hrxn.keys():
                data = []
                for ix in index:
                    if indiv[ix][db] is not None:
                        if rxn in odb.sset[lsset[index.index(ix)]].keys():
                            try:
                                data.append(indiv[ix][db][rxn][0])
                            except KeyError, e:
                                if failoninc:
                                    raise e
                                else:
                                    data.append(None)
                        else:
                            data.append(None)
                dbdat.append({'db': db,
                              'sys': str(rxn),
                              'color': odb.hrxn[rxn].color,
                              'data': data})
        mae = [errors[ix][self.dbse]['mae'] for ix in index]
        mape = [100 * errors[ix][self.dbse]['mape'] for ix in index]
        # form unique filename
        ixpre, ixsuf, ixmid = string_contrast(index)
        title = self.dbse + ' ' + ixpre + '[]' + ixsuf
        # generate matplotlib instructions and call or print
        try:
            import mpl
            import matplotlib.pyplot as plt
        except ImportError:
            # if not running from Canopy, print line to execute from Canopy
            print """filedict, htmlcode = mpl.threads(%s,\n    color='%s',\n    title='%s',\n    labels=%s,\n    mae=%s,\n    mape=%s\n    xlimit=%s\n    saveas=%s\n    mousetext=%s\n    mouselink=%s\n    mouseimag=%s\n    mousetitle=%s,\n    relpath=%s\n    graphicsformat=%s)\n\n""" % \
                (dbdat, color, title, ixmid, mae, mape, str(xlimit),
                repr(saveas), repr(mousetext), repr(mouselink), repr(mouseimag),
                repr(mousetitle), repr(relpath), repr(graphicsformat))
        else:
            # if running from Canopy, call mpl directly
            filedict, htmlcode = mpl.threads(dbdat, color=color, title=title, labels=ixmid, mae=mae, mape=mape, xlimit=xlimit, saveas=saveas, mousetext=mousetext, mouselink=mouselink, mouseimag=mouseimag, mousetitle=mousetitle, relpath=relpath, graphicsformat=graphicsformat)
            return filedict, htmlcode




    def table_generic(self, mtd, bas, columnplan, rowplan=['bas', 'mtd'],
        opt=['CP'], err=['mae'], sset=['tt'],
        benchmark='default', failoninc=True,
        landscape=False, standalone=True, subjoin=True,
        plotpath='', theme='', filename=None):
        """Prepares dictionary of errors for all combinations of *mtd*, *opt*,
        *bas* with respect to model chemistry *benchmark*, mindful of *failoninc*.
        Once error dictionary is ready, it and all other arguments are passed
        along to textables.table_generic.

        """
        # argument dbse=['DB4']  # TODO
        # gather list of model chemistries for table
        mcs = ['-'.join(prod) for prod in itertools.product(mtd, opt, bas)]

        # compute errors
        serrors = {}
        for mc in mcs:
            serrors[mc] = {}
            for ss in self.sset.keys():
                perr = self.compute_statistics(mc, benchmark=benchmark, sset=ss,
                    failoninc=failoninc, verbose=False, returnindiv=False)
                serrors[mc][ss] = {}
                serrors[mc][ss][self.dbse] = format_errors(perr[self.dbse], mode=3)
                for db in self.dbdict.keys():
                    serrors[mc][ss][db] = None if perr[db] is None else format_errors(perr[db], mode=3)

        textables.table_generic(dbse=dbse, serrors=serrors,
            mtd=mtd, bas=bas, columnplan=columnplan, rowplan=rowplan,
            opt=opt, err=err, sset=sset,
            landscape=landscape, standalone=standalone, subjoin=subjoin,
            plotpath=plotpath, theme=theme, filename=filename)

    def table_merge_abbr(self, mtd, bas, opt=['CP'], err=['mae'], benchmark='default', failoninc=True, plotpath='analysis/flats/mplflat_', theme='smmerge'):
        """Specialization of table_generic into table with minimal statistics
        (three S22 and three overall) plus embedded slat diagram as suitable
        for main paper. A single table is formed in sections by *bas* with
        lines *mtd* within each section.

        """
        rowplan = ['bas', 'mtd']
        columnplan = [
            ['l', r"""Method \& Basis Set""", '', textables.label, {}],
            ['d', r'S22', 'HB', textables.val, {'sset': 'hb', 'dbse': 'S22'}],
            ['d', r'S22', 'MX/DD', textables.val, {'sset': 'mxdd', 'dbse': 'S22'}],
            ['d', r'S22', 'TT', textables.val, {'sset': 'tt', 'dbse': 'S22'}],
            ['d', r'Overall', 'HB', textables.val, {'sset': 'hb', 'dbse': 'DB4'}],
            ['d', r'Overall', 'MX/DD', textables.val, {'sset': 'mxdd', 'dbse': 'DB4'}],
            ['d', r'Overall', 'TT', textables.val, {'sset': 'tt', 'dbse': 'DB4'}],
            ['l', r"""Error Distribution\footnotemark[1]""", r"""\includegraphics[width=6.67cm,height=3.5mm]{%s%s.pdf}""" % (plotpath, 'blank'), textables.graphics, {}],
            ['d', r'Time', '', textables.val, {'sset': 'tt-5min', 'dbse': 'NBC1'}]]
            # TODO Time column not right at all

        self.table_generic(mtd=mtd, bas=bas, columnplan=columnplan, rowplan=rowplan,
            opt=opt, err=err,
            benchmark=benchmark, failoninc=failoninc,
            landscape=False, standalone=True, subjoin=True,
            plotpath=plotpath, theme=theme, filename=None)
        # TODO: not handled: filename, TODO switch standalone

    def table_merge_suppmat(self, mtd, bas, opt=['CP'], err=['mae'], benchmark='default', failoninc=True, plotpath='analysis/flats/mplflat_', theme='lgmerge'):
        """Specialization of table_generic into table with as many statistics
        as will fit (mostly fullcurve and a few 5min) plus embedded slat
        diagram as suitable for supplementary material. Multiple tables are
        formed, one for each in *bas* with lines *mtd* within each table.

        """
        rowplan = ['bas', 'mtd']
        columnplan = [
            ['l', r"""Method \& Basis Set""", '', textables.label, {}],
            ['d', 'S22', 'HB', textables.val, {'sset': 'hb', 'dbse': 'S22'}],
            ['d', 'S22', 'MX', textables.val, {'sset': 'mx', 'dbse': 'S22'}],
            ['d', 'S22', 'DD', textables.val, {'sset': 'dd', 'dbse': 'S22'}],
            ['d', 'S22', 'TT', textables.val, {'sset': 'tt', 'dbse': 'S22'}],
            ['d', 'NBC10', 'MX', textables.val, {'sset': 'mx', 'dbse': 'NBC1'}],
            ['d', 'NBC10', 'DD', textables.val, {'sset': 'dd', 'dbse': 'NBC1'}],
            ['d', 'NBC10', 'TT', textables.val, {'sset': 'tt', 'dbse': 'NBC1'}],
            ['d', 'HBC6', 'HB/TT', textables.val, {'sset': 'tt', 'dbse': 'HBC1'}],
            ['d', 'HSG', 'HB', textables.val, {'sset': 'hb', 'dbse': 'HSG'}],
            ['d', 'HSG', 'MX', textables.val, {'sset': 'mx', 'dbse': 'HSG'}],
            ['d', 'HSG', 'DD', textables.val, {'sset': 'dd', 'dbse': 'HSG'}],
            ['d', 'HSG', 'TT', textables.val, {'sset': 'tt', 'dbse': 'HSG'}],
            ['d', 'Avg', 'TT ', textables.val, {'sset': 'tt', 'dbse': 'DB4'}],
            ['l', r"""Error Distribution\footnotemark[1]""", r"""\includegraphics[width=6.67cm,height=3.5mm]{%s%s.pdf}""" % (plotpath, 'blank'), textables.graphics, {}],
            ['d', 'NBC10', r"""TT\footnotemark[2]""", textables.val, {'sset': 'tt-5min', 'dbse': 'NBC1'}],
            ['d', 'HBC6', r"""TT\footnotemark[2] """, textables.val, {'sset': 'tt-5min', 'dbse': 'HBC1'}],
            ['d', 'Avg', r"""TT\footnotemark[2]""", textables.val, {'sset': 'tt-5min', 'dbse': 'DB4'}]]

        self.table_generic(mtd=mtd, bas=bas, columnplan=columnplan, rowplan=rowplan,
            opt=opt, err=err,
            benchmark=benchmark, failoninc=failoninc,
            landscape=True, standalone=True, subjoin=False,
            plotpath=plotpath, theme=theme, filename=None)
        # TODO: not handled: filename, TODO switch standalone


class DB4(Database):
    def __init__(self, pythonpath=None, loadfrompickle=False, path=None):
        """Initialize FourDatabases object from SuperDatabase"""
        Database.__init__(self, ['s22', 'nbc10', 'hbc6', 'hsg'], dbse='DB4', 
            pythonpath=pythonpath, loadfrompickle=loadfrompickle, path=path)

#        # load up data and definitions
#        self.load_qcdata_byproject('dft')
#        self.load_qcdata_byproject('pt2')
#        #self.load_qcdata_byproject('dhdft')
#        self.load_subsets()
        self.define_supersubsets()
        self.define_supermodelchems()

    def define_supersubsets(self):
        """

        """
        self.sset['tt'] = ['default', 'default', 'default', 'default']
        self.sset['hb'] = ['hb', None, 'default', 'hb']
        self.sset['mx'] = ['mx', 'mx', None, 'mx']
        self.sset['dd'] = ['dd', 'dd', None, 'dd']
        self.sset['mxdd'] = ['mxdd', 'default', None, 'mxdd']
        self.sset['pp'] = ['mxddpp', 'mxddpp', None, None]
        self.sset['np'] = ['mxddnp', 'mxddnp', None, 'mxdd']
        self.sset['tt-5min'] = ['default', '5min', '5min', 'default']
        self.sset['hb-5min'] = ['hb', None, '5min', 'hb']
        self.sset['mx-5min'] = ['mx', 'mx-5min', None, 'mx']
        self.sset['dd-5min'] = ['dd', 'dd-5min', None, 'dd']
        self.sset['mxdd-5min'] = ['mxdd', '5min', None, 'mxdd']
        self.sset['pp-5min'] = ['mxddpp', 'mxddpp-5min', None, None]
        self.sset['np-5min'] = ['mxddnp', 'mxddnp-5min', None, 'mxdd']

#    def benchmark(self):
#        """Returns the model chemistry label for the database's benchmark."""
#        return 'C2001BENCH'

    def define_supermodelchems(self):
        """

        """
        self.benchmark = 'C2011BENCH'
        self.mcs['C2011BENCH'] = ['S22A', 'NBC100', 'HBC60', 'HSG0']

        self.mcs['CCSD-CP-adz'] = ['CCSD-CP-adz', 'CCSD-CP-hadz', 'CCSD-CP-adz', 'CCSD-CP-hadz']
        self.mcs['CCSD-CP-atz'] = ['CCSD-CP-atz', 'CCSD-CP-hatz', 'CCSD-CP-atz', 'CCSD-CP-hatz']
        self.mcs['CCSD-CP-adtz'] = ['CCSD-CP-adtz', 'CCSD-CP-hadtz', 'CCSD-CP-adtz', 'CCSD-CP-hadtz']
        self.mcs['CCSD-CP-adtzadz'] = ['CCSD-CP-adtzadz', 'CCSD-CP-adtzhadz', 'CCSD-CP-adtzadz', 'CCSD-CP-adtzhadz']
        self.mcs['CCSD-CP-atzadz'] = ['CCSD-CP-atzadz', 'CCSD-CP-atzhadz', 'CCSD-CP-atzadz', 'CCSD-CP-atzhadz']
        self.mcs['CCSD-CP-atqzadz'] = ['CCSD-CP-atqzadz', 'CCSD-CP-atqzhadz', 'CCSD-CP-atqzadz', 'CCSD-CP-atqzhadz']
        self.mcs['CCSD-CP-atzadtz'] = ['CCSD-CP-atzadtz', 'CCSD-CP-atzhadtz', 'CCSD-CP-atzadtz', 'CCSD-CP-atzhadtz']
        self.mcs['CCSD-CP-atqzadtz'] = ['CCSD-CP-atqzadtz', 'CCSD-CP-atqzhadtz', 'CCSD-CP-atqzadtz', 'CCSD-CP-atqzhadtz']
        self.mcs['CCSD-CP-atqzatz'] = ['CCSD-CP-atqzatz', 'CCSD-CP-atqzhatz', 'CCSD-CP-atqzatz', 'CCSD-CP-atqzhatz']

        self.mcs['SCSCCSD-CP-adz'] = ['SCSCCSD-CP-adz', 'SCSCCSD-CP-hadz', 'SCSCCSD-CP-adz', 'SCSCCSD-CP-hadz']
        self.mcs['SCSCCSD-CP-atz'] = ['SCSCCSD-CP-atz', 'SCSCCSD-CP-hatz', 'SCSCCSD-CP-atz', 'SCSCCSD-CP-hatz']
        self.mcs['SCSCCSD-CP-adtz'] = ['SCSCCSD-CP-adtz', 'SCSCCSD-CP-hadtz', 'SCSCCSD-CP-adtz', 'SCSCCSD-CP-hadtz']
        self.mcs['SCSCCSD-CP-adtzadz'] = ['SCSCCSD-CP-adtzadz', 'SCSCCSD-CP-adtzhadz', 'SCSCCSD-CP-adtzadz', 'SCSCCSD-CP-adtzhadz']
        self.mcs['SCSCCSD-CP-atzadz'] = ['SCSCCSD-CP-atzadz', 'SCSCCSD-CP-atzhadz', 'SCSCCSD-CP-atzadz', 'SCSCCSD-CP-atzhadz']
        self.mcs['SCSCCSD-CP-atqzadz'] = ['SCSCCSD-CP-atqzadz', 'SCSCCSD-CP-atqzhadz', 'SCSCCSD-CP-atqzadz', 'SCSCCSD-CP-atqzhadz']
        self.mcs['SCSCCSD-CP-atzadtz'] = ['SCSCCSD-CP-atzadtz', 'SCSCCSD-CP-atzhadtz', 'SCSCCSD-CP-atzadtz', 'SCSCCSD-CP-atzhadtz']
        self.mcs['SCSCCSD-CP-atqzadtz'] = ['SCSCCSD-CP-atqzadtz', 'SCSCCSD-CP-atqzhadtz', 'SCSCCSD-CP-atqzadtz', 'SCSCCSD-CP-atqzhadtz']
        self.mcs['SCSCCSD-CP-atqzatz'] = ['SCSCCSD-CP-atqzatz', 'SCSCCSD-CP-atqzhatz', 'SCSCCSD-CP-atqzatz', 'SCSCCSD-CP-atqzhatz']

        self.mcs['SCSMICCSD-CP-adz'] = ['SCSMICCSD-CP-adz', 'SCSMICCSD-CP-hadz', 'SCSMICCSD-CP-adz', 'SCSMICCSD-CP-hadz']
        self.mcs['SCSMICCSD-CP-atz'] = ['SCSMICCSD-CP-atz', 'SCSMICCSD-CP-hatz', 'SCSMICCSD-CP-atz', 'SCSMICCSD-CP-hatz']
        self.mcs['SCSMICCSD-CP-adtz'] = ['SCSMICCSD-CP-adtz', 'SCSMICCSD-CP-hadtz', 'SCSMICCSD-CP-adtz', 'SCSMICCSD-CP-hadtz']
        self.mcs['SCSMICCSD-CP-adtzadz'] = ['SCSMICCSD-CP-adtzadz', 'SCSMICCSD-CP-adtzhadz', 'SCSMICCSD-CP-adtzadz', 'SCSMICCSD-CP-adtzhadz']
        self.mcs['SCSMICCSD-CP-atzadz'] = ['SCSMICCSD-CP-atzadz', 'SCSMICCSD-CP-atzhadz', 'SCSMICCSD-CP-atzadz', 'SCSMICCSD-CP-atzhadz']
        self.mcs['SCSMICCSD-CP-atqzadz'] = ['SCSMICCSD-CP-atqzadz', 'SCSMICCSD-CP-atqzhadz', 'SCSMICCSD-CP-atqzadz', 'SCSMICCSD-CP-atqzhadz']
        self.mcs['SCSMICCSD-CP-atzadtz'] = ['SCSMICCSD-CP-atzadtz', 'SCSMICCSD-CP-atzhadtz', 'SCSMICCSD-CP-atzadtz', 'SCSMICCSD-CP-atzhadtz']
        self.mcs['SCSMICCSD-CP-atqzadtz'] = ['SCSMICCSD-CP-atqzadtz', 'SCSMICCSD-CP-atqzhadtz', 'SCSMICCSD-CP-atqzadtz', 'SCSMICCSD-CP-atqzhadtz']
        self.mcs['SCSMICCSD-CP-atqzatz'] = ['SCSMICCSD-CP-atqzatz', 'SCSMICCSD-CP-atqzhatz', 'SCSMICCSD-CP-atqzatz', 'SCSMICCSD-CP-atqzhatz']

        self.mcs['CCSDT-CP-adz'] = ['CCSDT-CP-adz', 'CCSDT-CP-hadz', 'CCSDT-CP-adz', 'CCSDT-CP-hadz']
        self.mcs['CCSDT-CP-atz'] = ['CCSDT-CP-atz', 'CCSDT-CP-hatz', 'CCSDT-CP-atz', 'CCSDT-CP-hatz']
        self.mcs['CCSDT-CP-adtz'] = ['CCSDT-CP-adtz', 'CCSDT-CP-hadtz', 'CCSDT-CP-adtz', 'CCSDT-CP-hadtz']
        self.mcs['CCSDT-CP-adtzadz'] = ['CCSDT-CP-adtzadz', 'CCSDT-CP-adtzhadz', 'CCSDT-CP-adtzadz', 'CCSDT-CP-adtzhadz']
        self.mcs['CCSDT-CP-atzadz'] = ['CCSDT-CP-atzadz', 'CCSDT-CP-atzhadz', 'CCSDT-CP-atzadz', 'CCSDT-CP-atzhadz']
        self.mcs['CCSDT-CP-atqzadz'] = ['CCSDT-CP-atqzadz', 'CCSDT-CP-atqzhadz', 'CCSDT-CP-atqzadz', 'CCSDT-CP-atqzhadz']
        self.mcs['CCSDT-CP-atzadtz'] = ['CCSDT-CP-atzadtz', 'CCSDT-CP-atzhadtz', 'CCSDT-CP-atzadtz', 'CCSDT-CP-atzhadtz']
        self.mcs['CCSDT-CP-atqzadtz'] = ['CCSDT-CP-atqzadtz', 'CCSDT-CP-atqzhadtz', 'CCSDT-CP-atqzadtz', 'CCSDT-CP-atqzhadtz']
        self.mcs['CCSDT-CP-atqzatz'] = ['CCSDT-CP-atqzatz', 'CCSDT-CP-atqzhatz', 'CCSDT-CP-atqzatz', 'CCSDT-CP-atqzhatz']

    def make_pt2_flats(self):
        """Generate pieces for inclusion into tables for PT2 paper."""
        self.plot_all_flats(modelchem=None, sset='tt-5min', xlimit=4.0)

    def make_pt2_Figure_3(self):
        """Plot all the graphics needed for the calendar grey bars plot
        in Fig. 3 of PT2.

        """
        # Fig. bars (a)
        self.plot_bars(['MP2-CP-dz', 'MP2-CP-jadz', 'MP2-CP-hadz', 'MP2-CP-adz',
            'MP2-CP-tz', 'MP2-CP-matz', 'MP2-CP-jatz', 'MP2-CP-hatz', 'MP2-CP-atz',
            'MP2-CP-dtz', 'MP2-CP-jadtz', 'MP2-CP-hadtz', 'MP2-CP-adtz',
            'MP2-CP-qz', 'MP2-CP-aaqz', 'MP2-CP-maqz', 'MP2-CP-jaqz', 'MP2-CP-haqz', 'MP2-CP-aqz',
            'MP2-CP-tqz', 'MP2-CP-matqz', 'MP2-CP-jatqz', 'MP2-CP-hatqz', 'MP2-CP-atqz',
            'MP2-CP-a5z', 'MP2-CP-aq5z'])
        self.plot_bars(['SCSMP2-CP-dz', 'SCSMP2-CP-jadz', 'SCSMP2-CP-hadz', 'SCSMP2-CP-adz',
            'SCSMP2-CP-tz', 'SCSMP2-CP-matz', 'SCSMP2-CP-jatz', 'SCSMP2-CP-hatz', 'SCSMP2-CP-atz',
            'SCSMP2-CP-dtz', 'SCSMP2-CP-jadtz', 'SCSMP2-CP-hadtz', 'SCSMP2-CP-adtz',
            'SCSMP2-CP-qz', 'SCSMP2-CP-aaqz', 'SCSMP2-CP-maqz', 'SCSMP2-CP-jaqz', 'SCSMP2-CP-haqz', 'SCSMP2-CP-aqz',
            'SCSMP2-CP-tqz', 'SCSMP2-CP-matqz', 'SCSMP2-CP-jatqz', 'SCSMP2-CP-hatqz', 'SCSMP2-CP-atqz',
            'SCSMP2-CP-a5z', 'SCSMP2-CP-aq5z'])
        self.plot_bars(['SCSNMP2-CP-dz', 'SCSNMP2-CP-jadz', 'SCSNMP2-CP-hadz', 'SCSNMP2-CP-adz',
            'SCSNMP2-CP-tz', 'SCSNMP2-CP-matz', 'SCSNMP2-CP-jatz', 'SCSNMP2-CP-hatz', 'SCSNMP2-CP-atz',
            'SCSNMP2-CP-dtz', 'SCSNMP2-CP-jadtz', 'SCSNMP2-CP-hadtz', 'SCSNMP2-CP-adtz',
            'SCSNMP2-CP-qz', 'SCSNMP2-CP-aaqz', 'SCSNMP2-CP-maqz', 'SCSNMP2-CP-jaqz', 'SCSNMP2-CP-haqz', 'SCSNMP2-CP-aqz',
            'SCSNMP2-CP-tqz', 'SCSNMP2-CP-matqz', 'SCSNMP2-CP-jatqz', 'SCSNMP2-CP-hatqz', 'SCSNMP2-CP-atqz',
            'SCSNMP2-CP-a5z', 'SCSNMP2-CP-aq5z'])
        self.plot_bars([None, None, None, None,
            'SCSMIMP2-CP-tz', 'SCSMIMP2-CP-matz', 'SCSMIMP2-CP-jatz', 'SCSMIMP2-CP-hatz', 'SCSMIMP2-CP-atz',
            'SCSMIMP2-CP-dtz', 'SCSMIMP2-CP-jadtz', 'SCSMIMP2-CP-hadtz', 'SCSMIMP2-CP-adtz',
            'SCSMIMP2-CP-qz', 'SCSMIMP2-CP-aaqz', 'SCSMIMP2-CP-maqz', 'SCSMIMP2-CP-jaqz', 'SCSMIMP2-CP-haqz', 'SCSMIMP2-CP-aqz',
            'SCSMIMP2-CP-tqz', 'SCSMIMP2-CP-matqz', 'SCSMIMP2-CP-jatqz', 'SCSMIMP2-CP-hatqz', 'SCSMIMP2-CP-atqz',
            None, None])
        self.plot_bars(['DWMP2-CP-dz', 'DWMP2-CP-jadz', 'DWMP2-CP-hadz', 'DWMP2-CP-adz',
            'DWMP2-CP-tz', 'DWMP2-CP-matz', 'DWMP2-CP-jatz', 'DWMP2-CP-hatz', 'DWMP2-CP-atz',
            'DWMP2-CP-dtz', 'DWMP2-CP-jadtz', 'DWMP2-CP-hadtz', 'DWMP2-CP-adtz',
            'DWMP2-CP-qz', 'DWMP2-CP-aaqz', 'DWMP2-CP-maqz', 'DWMP2-CP-jaqz', 'DWMP2-CP-haqz', 'DWMP2-CP-aqz',
            'DWMP2-CP-tqz', 'DWMP2-CP-matqz', 'DWMP2-CP-jatqz', 'DWMP2-CP-hatqz', 'DWMP2-CP-atqz',
            'DWMP2-CP-a5z', 'DWMP2-CP-aq5z'])
        self.plot_bars(['MP2C-CP-dz', 'MP2C-CP-jadz', 'MP2C-CP-hadz', 'MP2C-CP-adz',
            'MP2C-CP-tz', 'MP2C-CP-matz', 'MP2C-CP-jatz', 'MP2C-CP-hatz', 'MP2C-CP-atz',
            'MP2C-CP-dtz', 'MP2C-CP-jadtz', 'MP2C-CP-hadtz', 'MP2C-CP-adtz',
            None, None, None, None, None, 'MP2C-CP-aqz',
            None, None, None, None, 'MP2C-CP-atqz',
            None, None])
        self.plot_bars(['MP2C-CP-atqzdz', 'MP2C-CP-atqzjadz', 'MP2C-CP-atqzhadz', 'MP2C-CP-atqzadz',
            'MP2C-CP-atqztz', 'MP2C-CP-atqzmatz', 'MP2C-CP-atqzjatz', 'MP2C-CP-atqzhatz', 'MP2C-CP-atqzatz',
            'MP2C-CP-atqzdtz', 'MP2C-CP-atqzjadtz', 'MP2C-CP-atqzhadtz', 'MP2C-CP-atqzadtz'])

        # Fig. bars (c)
        self.plot_bars(['MP2F12-CP-dz', 'MP2F12-CP-jadz', 'MP2F12-CP-hadz', 'MP2F12-CP-adz',
            'MP2F12-CP-tz', 'MP2F12-CP-matz', 'MP2F12-CP-jatz', 'MP2F12-CP-hatz', 'MP2F12-CP-atz',
            'MP2F12-CP-dtz', 'MP2F12-CP-jadtz', 'MP2F12-CP-hadtz', 'MP2F12-CP-adtz',
            'MP2F12-CP-aqz', 'MP2F12-CP-atqz'])
        self.plot_bars(['SCSMP2F12-CP-dz', 'SCSMP2F12-CP-jadz', 'SCSMP2F12-CP-hadz', 'SCSMP2F12-CP-adz',
            'SCSMP2F12-CP-tz', 'SCSMP2F12-CP-matz', 'SCSMP2F12-CP-jatz', 'SCSMP2F12-CP-hatz', 'SCSMP2F12-CP-atz',
            'SCSMP2F12-CP-dtz', 'SCSMP2F12-CP-jadtz', 'SCSMP2F12-CP-hadtz', 'SCSMP2F12-CP-adtz',
            'SCSMP2F12-CP-aqz', 'SCSMP2F12-CP-atqz'])
        self.plot_bars(['SCSNMP2F12-CP-dz', 'SCSNMP2F12-CP-jadz', 'SCSNMP2F12-CP-hadz', 'SCSNMP2F12-CP-adz',
            'SCSNMP2F12-CP-tz', 'SCSNMP2F12-CP-matz', 'SCSNMP2F12-CP-jatz', 'SCSNMP2F12-CP-hatz', 'SCSNMP2F12-CP-atz',
            'SCSNMP2F12-CP-dtz', 'SCSNMP2F12-CP-jadtz', 'SCSNMP2F12-CP-adtz', 'SCSNMP2F12-CP-adtz',
            'SCSNMP2F12-CP-aqz', 'SCSNMP2F12-CP-atqz'])
        self.plot_bars([None, None, None, None,
            'SCSMIMP2F12-CP-tz', 'SCSMIMP2F12-CP-matz', 'SCSMIMP2F12-CP-jatz', 'SCSMIMP2F12-CP-hatz', 'SCSMIMP2F12-CP-atz',
            'SCSMIMP2F12-CP-dtz', 'SCSMIMP2F12-CP-jadtz', 'SCSMIMP2F12-CP-hadtz', 'SCSMIMP2F12-CP-adtz',
            'SCSMIMP2F12-CP-aqz', 'SCSMIMP2F12-CP-atqz'])
        self.plot_bars(['DWMP2F12-CP-dz', 'DWMP2F12-CP-jadz', 'DWMP2F12-CP-hadz', 'DWMP2F12-CP-adz',
            'DWMP2F12-CP-tz', 'DWMP2F12-CP-matz', 'DWMP2F12-CP-jatz', 'DWMP2F12-CP-hatz', 'DWMP2F12-CP-atz',
            'DWMP2F12-CP-dtz', 'DWMP2F12-CP-jadtz', 'DWMP2F12-CP-hadtz', 'DWMP2F12-CP-adtz',
            'DWMP2F12-CP-aqz', 'DWMP2F12-CP-atqz'])
        self.plot_bars(['MP2CF12-CP-dz', 'MP2CF12-CP-jadz', 'MP2CF12-CP-hadz', 'MP2CF12-CP-adz',
            'MP2CF12-CP-tz', 'MP2CF12-CP-matz', 'MP2CF12-CP-jatz', 'MP2CF12-CP-hatz', 'MP2CF12-CP-atz',
            'MP2CF12-CP-dtz', 'MP2CF12-CP-jadtz', 'MP2CF12-CP-hadtz', 'MP2CF12-CP-adtz',
            'MP2CF12-CP-aqz', 'MP2CF12-CP-atqz'])
        self.plot_bars(['MP2CF12-CP-atqzdz', 'MP2CF12-CP-atqzjadz', 'MP2CF12-CP-atqzhadz', 'MP2CF12-CP-atqzadz',
            'MP2CF12-CP-atqztz', 'MP2CF12-CP-atqzmatz', 'MP2CF12-CP-atqzjatz', 'MP2CF12-CP-atqzhatz', 'MP2CF12-CP-atqzatz',
            'MP2CF12-CP-atqzdtz', 'MP2CF12-CP-atqzjadtz', 'MP2CF12-CP-atqzhadtz', 'MP2CF12-CP-atqzadtz'])

    def make_pt2_Figure_2(self):
        """Plot all the graphics needed for the diffuse augmented grey
        bars plot in Fig. 2 of PT2.

        """
        # Fig. bars (a)
        self.plot_bars(['MP2-CP-adz', 'MP2-CP-atz', 'MP2-CP-adtz',
            'MP2-CP-aqz', 'MP2-CP-atqz', 'MP2-CP-a5z', 'MP2-CP-aq5z'])
        self.plot_bars(['SCSMP2-CP-adz', 'SCSMP2-CP-atz',
            'SCSMP2-CP-adtz', 'SCSMP2-CP-aqz', 'SCSMP2-CP-atqz',
            'SCSMP2-CP-a5z', 'SCSMP2-CP-aq5z'])
        self.plot_bars(['SCSNMP2-CP-adz', 'SCSNMP2-CP-atz',
            'SCSNMP2-CP-adtz', 'SCSNMP2-CP-aqz', 'SCSNMP2-CP-atqz',
            'SCSNMP2-CP-a5z', 'SCSNMP2-CP-aq5z'])
        self.plot_bars(['SCSMIMP2-CP-atz', 'SCSMIMP2-CP-atz',
            'SCSMIMP2-CP-adtz', 'SCSMIMP2-CP-aqz', 'SCSMIMP2-CP-atqz'])
        self.plot_bars(['SCSMIMP2-CP-tz', 'SCSMIMP2-CP-tz',
            'SCSMIMP2-CP-dtz', 'SCSMIMP2-CP-qz', 'SCSMIMP2-CP-tqz'])
        self.plot_bars(['DWMP2-CP-adz', 'DWMP2-CP-atz', 'DWMP2-CP-adtz',
            'DWMP2-CP-aqz', 'DWMP2-CP-atqz', 'DWMP2-CP-a5z', 'DWMP2-CP-aq5z'])
        self.plot_bars(['MP2C-CP-adz', 'MP2C-CP-adtzadz',
            'MP2C-CP-atqzadz', 'MP2C-CP-aq5zadz', 'MP2C-CP-atz',
            'MP2C-CP-atqzatz', 'MP2C-CP-aq5zatz', 'MP2C-CP-adtz',
            'MP2C-CP-atqzadtz', 'MP2C-CP-aqz', 'MP2C-CP-atqz'])

        # Fig. bars (b)
        self.plot_bars(['MP3-CP-adz', 'MP3-CP-adtzadz', 'MP3-CP-atqzadz',
            'MP3-CP-atz', 'MP3-CP-atqzatz', 'MP3-CP-adtz', 'MP3-CP-atqzadtz'])
        self.plot_bars(['MP25-CP-adz', 'MP25-CP-adtzadz', 'MP25-CP-atqzadz',
            'MP25-CP-atz', 'MP25-CP-atqzatz', 'MP25-CP-adtz', 'MP25-CP-atqzadtz'])
        self.plot_bars(['CCSD-CP-adz', 'CCSD-CP-adtzadz', 'CCSD-CP-atqzadz',
            'CCSD-CP-atz', 'CCSD-CP-atqzatz', 'CCSD-CP-adtz', 'CCSD-CP-atqzadtz'])
        self.plot_bars(['SCSCCSD-CP-adz', 'SCSCCSD-CP-adtzadz',
            'SCSCCSD-CP-atqzadz', 'SCSCCSD-CP-atz', 'SCSCCSD-CP-atqzatz',
            'SCSCCSD-CP-adtz', 'SCSCCSD-CP-atqzadtz'])
        self.plot_bars(['SCSMICCSD-CP-adz', 'SCSMICCSD-CP-adtzadz',
            'SCSMICCSD-CP-atqzadz', 'SCSMICCSD-CP-atz', 'SCSMICCSD-CP-atqzatz',
            'SCSMICCSD-CP-adtz', 'SCSMICCSD-CP-atqzadtz'])
        self.plot_bars(['CCSDT-CP-adz', 'CCSDT-CP-adtzadz',
            'CCSDT-CP-atqzadz', 'CCSDT-CP-atz', 'CCSDT-CP-atqzatz',
            'CCSDT-CP-adtz', 'CCSDT-CP-atqzadtz'])

        # Fig. bars (c)
        self.plot_bars(['MP2F12-CP-adz', 'MP2F12-CP-atz', 'MP2F12-CP-adtz',
            'MP2F12-CP-aqz', 'MP2F12-CP-atqz'])
        self.plot_bars(['SCSMP2F12-CP-adz', 'SCSMP2F12-CP-atz',
            'SCSMP2F12-CP-adtz', 'SCSMP2F12-CP-aqz', 'SCSMP2F12-CP-atqz'])
        self.plot_bars(['SCSNMP2F12-CP-adz', 'SCSNMP2F12-CP-atz',
            'SCSNMP2F12-CP-adtz', 'SCSNMP2F12-CP-aqz',
            'SCSNMP2F12-CP-atqz'])
        self.plot_bars(['SCSMIMP2F12-CP-atz', 'SCSMIMP2F12-CP-atz',
            'SCSMIMP2F12-CP-adtz', 'SCSMIMP2F12-CP-aqz',
            'SCSMIMP2F12-CP-atqz'])
        self.plot_bars(['SCSMIMP2F12-CP-tz', 'SCSMIMP2F12-CP-tz', 'SCSMIMP2F12-CP-dtz'])
        self.plot_bars(['DWMP2F12-CP-adz', 'DWMP2F12-CP-atz',
            'DWMP2F12-CP-adtz', 'DWMP2F12-CP-aqz', 'DWMP2F12-CP-atqz'])
        self.plot_bars(['MP2CF12-CP-adz', 'MP2CF12-CP-adtzadz',
            'MP2CF12-CP-atqzadz', 'MP2CF12-CP-atz', 'MP2CF12-CP-atqzatz',
            'MP2CF12-CP-adtz', 'MP2CF12-CP-atqzadtz', 'MP2CF12-CP-aqz',
            'MP2CF12-CP-atqz'])

        # Fig. bars (d)
        self.plot_bars(['CCSDAF12-CP-adz', 'CCSDAF12-CP-adtzadz', 'CCSDAF12-CP-atqzadz'])
        self.plot_bars(['CCSDBF12-CP-adz', 'CCSDBF12-CP-adtzadz', 'CCSDBF12-CP-atqzadz'])
        self.plot_bars(['SCSCCSDAF12-CP-adz', 'SCSCCSDAF12-CP-adtzadz', 'SCSCCSDAF12-CP-atqzadz'])
        self.plot_bars(['SCSCCSDBF12-CP-adz', 'SCSCCSDBF12-CP-adtzadz', 'SCSCCSDBF12-CP-atqzadz'])
        self.plot_bars(['SCMICCSDAF12-CP-adz', 'SCMICCSDAF12-CP-adtzadz', 'SCMICCSDAF12-CP-atqzadz'])
        self.plot_bars(['SCMICCSDBF12-CP-adz', 'SCMICCSDBF12-CP-adtzadz', 'SCMICCSDBF12-CP-atqzadz'])
        self.plot_bars(['CCSDTAF12-CP-adz', 'CCSDTAF12-CP-adtzadz', 'CCSDTAF12-CP-atqzadz'])
        self.plot_bars(['CCSDTBF12-CP-adz', 'CCSDTBF12-CP-adtzadz', 'CCSDTBF12-CP-atqzadz'])
        self.plot_bars(['DWCCSDTF12-CP-adz', 'DWCCSDTF12-CP-adtzadz', 'DWCCSDTF12-CP-atqzadz'])


class ThreeDatabases(Database):
    """

    """
    def __init__(self, pythonpath=None):
        """Initialize ThreeDatabases object from Database"""
        Database.__init__(self, ['s22', 'a24', 'hsg'], dbse='DB3', pythonpath=None)

        # load up data and definitions
        self.load_qcdata_byproject('pt2')
        self.load_qcdata_byproject('dilabio')
        self.load_qcdata_byproject('f12dilabio')
        self.load_subsets()
        self.define_supersubsets()
        self.define_supermodelchems()

    def define_supersubsets(self):
        """

        """
        self.sset['tt'] = ['default', 'default', 'default']
        self.sset['hb'] = ['hb', 'hb', 'hb']
        self.sset['mx'] = ['mx', 'mx', 'mx']
        self.sset['dd'] = ['dd', 'dd', 'dd']
        self.sset['mxdd'] = ['mxdd', 'mxdd', 'mxdd']
        self.sset['pp'] = ['mxddpp', 'mxddpp', 'mxddpp']
        self.sset['np'] = ['mxddnp', 'mxddnp', 'mxddnp']
        self.sset['tt-5min'] = ['default', 'default', 'default']
        self.sset['hb-5min'] = ['hb', 'hb', 'hb']
        self.sset['mx-5min'] = ['mx', 'mx', 'mx']
        self.sset['dd-5min'] = ['dd', 'dd', 'dd']
        self.sset['mxdd-5min'] = ['mxdd', 'mxdd', 'mxdd']
        self.sset['pp-5min'] = ['mxddpp', 'mxddpp', 'mxddpp']
        self.sset['np-5min'] = ['mxddnp', 'mxddnp', 'mxddnp']
        self.sset['weak'] = ['weak', 'weak', 'weak']
        self.sset['weak_hb'] = ['weak_hb', None, 'weak_hb']
        self.sset['weak_mx'] = ['weak_mx', 'weak_mx', 'weak_mx']
        self.sset['weak_dd'] = ['weak_dd', 'weak_dd', 'weak_dd']

    def define_supermodelchems(self):
        """

        """
        self.mc['CCSD-CP-adz'] = ['CCSD-CP-adz', 'CCSD-CP-hadz', 'CCSD-CP-adz']
        self.mc['CCSD-CP-atz'] = ['CCSD-CP-atz', 'CCSD-CP-hatz', 'CCSD-CP-atz']
        self.mc['CCSD-CP-adtz'] = ['CCSD-CP-adtz', 'CCSD-CP-hadtz', 'CCSD-CP-adtz']
        self.mc['CCSD-CP-adtzadz'] = ['CCSD-CP-adtzadz', 'CCSD-CP-adtzhadz', 'CCSD-CP-adtzadz']
        self.mc['CCSD-CP-atzadz'] = ['CCSD-CP-atzadz', 'CCSD-CP-atzhadz', 'CCSD-CP-atzadz']
        self.mc['CCSD-CP-atqzadz'] = ['CCSD-CP-atqzadz', 'CCSD-CP-atqzhadz', 'CCSD-CP-atqzadz']
        self.mc['CCSD-CP-atzadtz'] = ['CCSD-CP-atzadtz', 'CCSD-CP-atzhadtz', 'CCSD-CP-atzadtz']
        self.mc['CCSD-CP-atqzadtz'] = ['CCSD-CP-atqzadtz', 'CCSD-CP-atqzhadtz', 'CCSD-CP-atqzadtz']
        self.mc['CCSD-CP-atqzatz'] = ['CCSD-CP-atqzatz', 'CCSD-CP-atqzhatz', 'CCSD-CP-atqzatz']

        self.mc['CCSDT-CP-adz'] = ['CCSDT-CP-adz', 'CCSDT-CP-hadz', 'CCSDT-CP-adz']
        self.mc['CCSDT-CP-atz'] = ['CCSDT-CP-atz', 'CCSDT-CP-hatz', 'CCSDT-CP-atz']
        self.mc['CCSDT-CP-adtz'] = ['CCSDT-CP-adtz', 'CCSDT-CP-hadtz', 'CCSDT-CP-adtz']
        self.mc['CCSDT-CP-adtzadz'] = ['CCSDT-CP-adtzadz', 'CCSDT-CP-adtzhadz', 'CCSDT-CP-adtzadz']
        self.mc['CCSDT-CP-atzadz'] = ['CCSDT-CP-atzadz', 'CCSDT-CP-atzhadz', 'CCSDT-CP-atzadz']
        self.mc['CCSDT-CP-atqzadz'] = ['CCSDT-CP-atqzadz', 'CCSDT-CP-atqzhadz', 'CCSDT-CP-atqzadz']
        self.mc['CCSDT-CP-atzadtz'] = ['CCSDT-CP-atzadtz', 'CCSDT-CP-atzhadtz', 'CCSDT-CP-atzadtz']
        self.mc['CCSDT-CP-atqzadtz'] = ['CCSDT-CP-atqzadtz', 'CCSDT-CP-atqzhadtz', 'CCSDT-CP-atqzadtz']
        self.mc['CCSDT-CP-atqzatz'] = ['CCSDT-CP-atqzatz', 'CCSDT-CP-atqzhatz', 'CCSDT-CP-atqzatz']

# print certain statistic for all 4 db and summary and indiv sys if min or max
